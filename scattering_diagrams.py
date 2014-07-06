from sage.misc.cachefunc import cached_method                                                                                                                             

class ScatteringWall(SageObject):
    """r
        A class to encode a wall in a scattering diagram
    """
    def __init__(self, slope, f):
        self.slope = tuple(slope)
        self.function = f
    
    def __hash__(self):
        return hash((self.slope, self.function))
    
    def __repr__(self):
        return "A wall of slope " + str(self.slope) + " with scattering function " + str(self.function)

    def _latex_(self):
        pass

    def tikz(self, radius=10,show_label=False):
        output = "\\draw[color=black,line width=1pt] (0,0) -- "
        endpoint = -N(radius * vector(self.slope)/norm(vector(self.slope)))
        output += str(endpoint) + ";\n"
        # Do something here to show funtions
        return output

class ScatteringDiagram(SageObject):
    """r
        A class to hold the scattering diagram (usually up to a certain depth)
    """
    def __init__(self, walls):
        self.walls = tuple(walls)

    def __hash__(self):
        return hash(self.walls)
    
    def __repr__(self):
        return "A scattering diagram with " + str(len(self.walls)) + " walls in it."

    def __getitem__(self, idx):
        return self.walls[idx]

    def __setitem__(self, idx, wall):
        new_diagram = list(self.walls)
        new_diagram[idx] = wall
        self.walls = tuple(new_diagram)

    def __len__(self):
        return len(self.walls)

    def _latex_(self):
        pass

    def insert(self, position, wall):
        new_diagram = list(self.walls)
        new_diagram.insert(position, wall)
        self.walls = tuple(new_diagram)

    def append(self, wall):
        new_diagram = list(self.walls)
        new_diagram.append(wall)
        self.walls = tuple(new_diagram)

    def tikz(self, radius=10,show_label=False):
        output = ""
        for W in self.walls:
            output += W.tikz(radius=radius, show_label=show_label)
        return output

class ScatteringRing(SageObject):
    """r
        A rank-2 scattering ring
    """
    def __init__(self, b, c):
        self._b = b
        self._c = c
        self.scatter_ring = LaurentPolynomialRing(QQ, 'x,y')
        self.scatter_field = self.scatter_ring.fraction_field()
        (self.x, self.y) = self.scatter_ring.gens()
    
    @cached_method
    def scattering_diagram(self, max_degree):
        # this map should be made recursive
        W = ScatteringWall
        x = self.x
        y = self.y
        b = self._b 
        c = self._c
        diagram = ScatteringDiagram([ W((0,-1),1+y**c), W((1,0),1+x**b), W((0,1),1+y**c), W((-1,0),1+x**b) ])
        for k in range(2,max_degree+1):
           diagram = self._add_walls(diagram, k)
        return diagram

    def _path_product(self, scatter_diagram, k):
        X = self.scatter_field(self.x)
        for wall in scatter_diagram:
            slope = wall.slope
            f = wall.function
            g = gcd(slope)
            new_x = self.scatter_field(self.x*f**(-slope[1]/g))
            new_y = self.scatter_field(self.y*f**(slope[0]/g))
            X = X(x=new_x,y=new_y)
            X = taylor(X,(self.x,0),(self.y,0),k-1)
        return self.scatter_ring(X)

    def _add_walls(self, diagram, k):
        P = self._path_product(diagram, k)
        for i in range(2,k-1):
            g = gcd(i-1, k-i-1)
            C = g*P.coefficient(self.x**i*self.y**(k-1-i))/(k-1-i)
            if C != 0:
                for j in range(2,len(diagram)):
                    wall = diagram[j]
                    slope = wall.slope
                    f = wall.function
                    if slope[0] == (i-1)/g and slope[1] == (k-1-i)/g:
                        replacement_polynomial = f+C*self.x**(i-1)*self.y**(k-1-i)
                        #replacement_polynomial = replacement_polynomial.polynomial(self.scatter_ring)
                        replacement_wall = ScatteringWall( slope, replacement_polynomial )
                        diagram[j] = replacement_wall 
                        break
                    elif slope[0] == 0 or (i-1)*slope[1] > (k-1-i)*slope[0]*sign(slope[0]):
                        new_polynomial = 1+C*self.x**(i-1)*self.y**(k-1-i)
                        #new_polynomial = new_polynomial.polynomial(self.scatter_ring)
                        new_wall = ScatteringWall( ((i-1)/g,(k-1-i)/g), new_polynomial )
                        diagram.insert(j, new_wall)
                        break
                    else:
                        j += 1
        return diagram

    @cached_method
    def broken_lines(self, init_momentum, end_point):
        degree = _monomial_degree(init_momentum)
        if any( x > 0 for x in degree ):
            return (BrokenLine(init_momentum),)

        total_degree = -sum(degree)
        diagram = self.scattering_diagram(total_degree)
        
        clockwise_diagram = ScatteringDiagram([ W for W in diagram if _side(W.slope,degree) == "clockwise"  ])
        counterclockwise_diagram = ScatteringDiagram(reversed([ W for W in diagram if _side(W.slope,degree) == "counterclockwise"  ]))
        print clockwise_diagram.walls
        print counterclockwise_diagram.walls

        clockwise_lines = [BrokenLine(init_momentum)]
        for wall in clockwise_diagram:
            temp_lines = []
            for line in clockwise_lines:
                temp_lines.append(line)
                last_momentum = line.monomial()
                for new_segment in self._scatter_at_wall(last_momentum,wall):
                    if _monomial_degree(new_segment.monomial)[1]<0 and new_segment.monomial != last_momentum:
                        new_line=copy(line)
                        new_line.append_segment(new_segment)
                        temp_lines.append(new_line)
            clockwise_lines = temp_lines

        counterclockwise_lines = [BrokenLine(init_momentum)]
        for wall in counterclockwise_diagram:
            temp_lines = []
            for line in counterclockwise_lines:
                temp_lines.append(line)
                last_momentum = line.monomial()
                for new_segment in self._scatter_at_wall(last_momentum,wall):
                    if _monomial_degree(new_segment.monomial)[0]<0 and new_segment.monomial != last_momentum:
                        new_line=copy(line)
                        new_line.append_segment(new_segment)
                        temp_lines.append(new_line)
            counterclockwise_lines = temp_lines
        
        broken_lines = []
        for line in clockwise_lines:
            final_direction = _monomial_degree(line.monomial())
            if final_direction[0] > 0:
                broken_lines.append(line)
            elif end_point[0]*final_direction[1]-end_point[1]*final_direction[0] < 0:
                broken_lines.append(line)
        for line in counterclockwise_lines:
            final_direction = _monomial_degree(line.monomial())
            if final_direction[1] > 0:
                broken_lines.append(line)
            elif end_point[0]*final_direction[1]-end_point[1]*final_direction[0] > 0:
                broken_lines.append(line)
        
        return tuple(broken_lines)
    
    def _scatter_at_wall(self, momentum, wall):
        momentum_exp = _monomial_degree(momentum)
        exp = abs(wall.slope[0]*momentum_exp[1]-wall.slope[1]*momentum_exp[0])
        scattered_function = momentum*(wall.function**exp)
        scattered_monomials = scattered_function.monomials()
        scattered_segments = []
        for mon in scattered_monomials:
            scattered_segments.append(BrokenLineSegment(scattered_function.monomial_coefficient(mon), mon, wall.slope))
        return tuple(scattered_segments)

class BrokenLine(SageObject):
    """r
        A broken line
    """
    
    def __init__(self, initial_momentum):
        if type(initial_momentum) == tuple:
            # do something to pass integers instead of monomials
            pass
        self.line_segments = (BrokenLineSegment(1, initial_momentum, None),)

    def __repr__(self):
        output = "A broken line with "
        output += str(len(self.line_segments))
        output += " segment"
        if len(self.line_segments) > 1:
            output +="s"
            output +=", initial momentum "
            output += str(self.line_segments[0].monomial)
            output += ", final "
        else:
            output += ", "
        output += "momentum "
        output += str(self.line_segments[-1].monomial)
        output += ", and coefficient "
        output += str(self.coefficient())
        output +="."
        return output

    def append_segment(self, segment):
        self.line_segments += (segment,)

    def coefficient(self):
        return prod([ x.coefficient for x in self.line_segments ])

    def monomial(self):
        return self.line_segments[-1].monomial

    def tikz(self, end_point, radius):
        output = ""
        current_point = end_point
        for segment in reversed(self.line_segments):
            current_direction  = _monomial_degree(segment.monomial)
            intersection_slope = segment.scattering_wall
            if intersection_slope != None:
                var('s','t')
                solution_dict = solve([current_direction[1]*(s-current_point[0])==current_direction[0]*(t-current_point[1]),intersection_slope[1]*s==intersection_slope[0]*t],(s,t),solution_dict=True)[0]
                final_point = (solution_dict[s],solution_dict[t])
            else:
                final_point = (current_point[0]+radius*current_direction[0],current_point[1]+radius*current_direction[1])
            output += "  \\draw[color=red,line width=1pt] ("
            output += str(current_point[0])+","+str(current_point[1])+") -- ("
            output += str(final_point[0])+","+str(final_point[1])+");\n"
            current_point = final_point
            if intersection_slope != None:
                output += "  \\draw[color=red,fill=red] ("+str(final_point[0])+","+str(final_point[1])+") circle (3pt);\n"
        return output


class BrokenLineSegment(SageObject):
    """r
        A segment of broken line
    """

    def __init__(self, coefficient, monomial, scattering_wall):
        self.coefficient = coefficient
        self.monomial = monomial
        self.scattering_wall = scattering_wall

    def __repr__(self):
        output = "A broken line segment with coefficient " 
        output += str(self.coefficient)
        output += ", slope "
        output += str(self.monomial)
        if self.scattering_wall != None:
            output += ", scattering on the wall of slope "
            output += str(self.scattering_wall)
        output += "."
        return output

###
# helper functions
##

def _monomial_degree(monomial):
        monomial = monomial.parent().fraction_field()(monomial)
        return vector(monomial.numerator().exponents()[0])-vector(monomial.denominator().exponents()[0])

def _side(to_check, reference):
    if to_check == (-1,0) or to_check == (0,1):
        return "clockwise"
    if  to_check == (0,-1) or to_check == (1,0):
        return "counterclockwise"
    if reference[0]*to_check[1] < reference[1]*to_check[0]:
        return "clockwise"
    if reference[0]*to_check[1] > reference[1]*to_check[0]:
        return "counterclockwise"
    return "colinear"

#################################
# Old stuff
################################


#    #scattering methods
#        
#        
#
#    def find_walls(self,init_momentum):
#        counterclockwise_root_walls = [(0,1)]
#        clockwise_root_walls = [(1,0)]
#        counterclockwise_split_root_walls = []
#        clockwise_split_root_walls = []
#        
#        momentum_exp = vector(init_momentum.numerator().exponents()[0])-vector(init_momentum.denominator().exponents()[0])
#        check1 = true
#        n=0
#        while check1:
#            root = (sqrt((sqrt(self._c/self._b)**(n%2)*chebyshev_U(n,sqrt(self._b*self._c)/2))**2),sqrt((sqrt(self._b/self._c)**((n-1)%2)*chebyshev_U(n-1,sqrt(self._b*self._c)/2))**2))
#            check1 = false
#            if momentum_exp[0]*root[1]-momentum_exp[1]*root[0]>0:
#                if momentum_exp[0]+self._b*root[0]<0:
#                    counterclockwise_root_walls = [root]+counterclockwise_root_walls
#                    check1 = true
#            elif momentum_exp[0]*root[1]-momentum_exp[1]*root[0]<0:
#                if momentum_exp[1]+self._c*root[1]<0:
#                    counterclockwise_split_root_walls.append(root)
#                    check1 = true
#            n += 1
#        check2 = true
#        n=0
#        while check2:
#            root = (sqrt((sqrt(self._c/self._b)**((n-1)%2)*chebyshev_U(n-1,sqrt(self._b*self._c)/2))**2),sqrt((sqrt(self._b/self._c)**(n%2)*chebyshev_U(n,sqrt(self._b*self._c)/2))**2))
#            check2 = false
#            if momentum_exp[0]*root[1]-momentum_exp[1]*root[0]<0:
#                if momentum_exp[1]+self._c*root[1]<0:
#                    clockwise_root_walls = [root]+clockwise_root_walls
#                    check2 = true
#            elif momentum_exp[0]*root[1]-momentum_exp[1]*root[0]>0:
#                if momentum_exp[0]+self._b*root[0]<0:
#                    clockwise_split_root_walls.append(root)
#                    check2 = true
#            n += 1
#
#        ###########
#        # Imaginary walls needed
#        ###########
#        left_imaginary_walls = []
#        right_imaginary_walls = []
#        counterclockwise_walls = clockwise_split_root_walls+left_imaginary_walls+counterclockwise_root_walls
#        clockwise_walls = counterclockwise_split_root_walls+right_imaginary_walls+clockwise_root_walls
#
#        #print "clockwise=",clockwise_walls
#        #print "counterclockwise_walls=",counterclockwise_walls
#
#        return [clockwise_walls,counterclockwise_walls]
#        
#        
#
#
#    def GreedyElementBrokenLine(self,a1,a2,Q):
#        """
#        Output:
#            -Returns the greedy basis element
#        """
#
#        output = 0
#        for line in self.broken_lines(self._x**(-a1)*self._y**(-a2),Q):
#            coeff = prod([i[0] for i in line])
#            output += coeff*line[-1][1]
#        return output 
#        
#    
#    def draw_broken_line(self,line,Q,rgb=(1,0,0)):
#        color = "{rgb:red,"+str(rgb[0])+";green,"+str(rgb[1])+";blue,"+str(rgb[2])+"}"
#        output = ""
#        current_point = Q
#        for i in range(1,line.__len__()):
#            current_direction  = self.monomial_degree(line[-i][1])
#            intersection_slope = line[-i][2]
#            #print "line=",line[-i]
#            #print "direction=",current_direction
#            #print "intersection=",intersection_slope
#            var('s','t')
#            solution_dict = solve([current_direction[1]*(s-current_point[0])==current_direction[0]*(t-current_point[1]),intersection_slope[1]*s==intersection_slope[0]*t],(s,t),solution_dict=True)[0]
#            final_point = (solution_dict[s],solution_dict[t])
#            output += "  \\draw[color="+color+",line width=1pt] ("
#            output += str(current_point[0])+","+str(current_point[1])+") -- ("
#            output += str(final_point[0])+","+str(final_point[1])+");\n"
#            output += "  \\draw[color="+color+",fill="+color+"] ("+str(final_point[0])+","+str(final_point[1])+") circle (3pt);\n"
#            current_point = final_point
#        current_direction  = self.monomial_degree(line[0][1])
#        final_point = (current_point[0]+max(Q)*current_direction[0],current_point[1]+max(Q)*current_direction[1])
#        output += "  \\draw[color="+color+",line width=1pt] ("
#        output += str(current_point[0])+","+str(current_point[1])+") -- ("
#        output += str(final_point[0])+","+str(final_point[1])+");\n"
#        return output
#
#
#    def draw_walls(self,a1,a2,Q):
#        output = ""
#        walls = self.find_walls(self._x**(-a1)*self._y**(-a2))
#        length = ceil(sqrt(Q[0]**2+Q[1]**2))
#        for w in walls[0][:-1]+walls[1][:-1]:
#            wall_length = ceil(sqrt(w[0]**2+w[1]**2))
#            output += "  \\draw[color=black,line width=1pt] (0,0) -- ("
#            output += str(-ceil(self._b*length*w[0]/wall_length))+","+str(-ceil(self._c*length*w[1]/wall_length))+");\n"
#        for w in [(0,1),(1,0)]:
#            output += "  \\draw[color=black,line width=1pt] (0,0) -- ("
#            output += str(self._b*length*w[0])+","+str(self._c*length*w[1])+");\n"
#        #print output
#        return output
#        
#        
#    def draw_broken_lines(self,a1,a2,Q):
#        #creates tikz pictures of all broken lines with initial momentum x**{-a1}y**{-a2} ending at Q
#        working_dir="/tmp/"
#        #filename = "broken_lines"+str(self._b)+str(self._c)+"-("+str(a1)+","+str(a2)+")"#-("+str(Q[0])+","+str(Q[1])+")"
#        filename = "sage_output"
#        filename += ".tex"
#        print filename
#        
#        TeXFile=open(working_dir+filename,'w')
#        TeXFile.write("\\documentclass{article}\n")
#        TeXFile.write("\\usepackage{amsmath, amssymb, latexsym, tikz}\n")
#        TeXFile.write("\\usepackage{pgflibraryarrows,pgflibrarysnakes}\n\n")
#        TeXFile.write("\\begin{document}\n\n")
#        TeXFile.write("Let $b="+str(self._b)+"$ and $c="+str(self._c)+"$.  ")
#        TeXFile.write("We consider broken lines with initial momentum $")
#        if a1 != 0:
#            TeXFile.write("x")
#            if a1 != -1:
#                TeXFile.write("**{")
#                if a1 > 0:
#                    TeXFile.write("-")
#                    TeXFile.write(str(a1)+"}")
#                else:
#                    TeXFile.write(str(-a1)+"}")
#        if a2 != 0:
#            TeXFile.write("y")
#            if a2 != -1:
#                TeXFile.write("**{")
#                if a2 > 0:
#                    TeXFile.write("-")
#                    TeXFile.write(str(a2)+"}")
#                else:
#                    TeXFile.write(str(-a2)+"}")
#        
#        TeXFile.write("$ which end at $Q=("+str(Q[0])+","+str(Q[1])+")$.  These are given by:\\\\\n\n")
#        TeXFile.write("\\resizebox{4.5in}{4.5in}{\n")
#        TeXFile.write(" \\begin{tikzpicture}\n")
#        TeXFile.write(self.draw_walls(a1,a2,Q))
#        
#        for line in self.broken_lines(self._x**(-a1)*self._y**(-a2),Q):
#            TeXFile.write(self.draw_broken_line(line,Q))
#
#
#        TeXFile.write("  \\draw[fill=black] ("+str(Q[0])+","+str(Q[1])+") circle (3pt);\n")
#        
#        TeXFile.write(" \\end{tikzpicture}}\\\\\n\n")
#            
#        TeXFile.write("\\end{document}")
#        TeXFile.close()
#        
#        import subprocess
#        subprocess.call(['pdflatex', '-halt-on-error', filename], cwd=working_dir, stdout=subprocess.PIPE)
