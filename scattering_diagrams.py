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

    def show(self, radius=10,show_label=False):
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

    def show(self, radius=10,show_label=False):
        output = ""
        for W in self.walls:
            output += W.show(radius=radius, show_label=show_label)
        return output

class ScatteringRing(SageObject):
    """r
        A rank-2 scattering ring
    """
    def __init__(self, b, c):
        self._b = b
        self._c = c
        
        var('x,y')
        self.scatter_ring=ZZ[x,y].fraction_field()
        self.x = self.scatter_ring(x)
        self.y = self.scatter_ring(y)
    
    @cached_method
    def scattering_diagram(self, depth):
        # this map should be made recursive
        W = ScatteringWall
        x = self.x
        y = self.y
        b = self._b 
        c = self._c
        diagram = ScatteringDiagram([ W((-1,0),1+x**b), W((0,-1),1+y**c), W((1,0),1+x**b), W((0,1),1+y**c) ])
        for k in range(2,depth+1):
           diagram = self._add_walls(diagram, k)
        return diagram

    def _path_product(self, scatter_diagram, k):
        X = self.x
        for wall in scatter_diagram:
            slope = wall.slope
            f = wall.function
            g = gcd(slope)
            X = X(x=self.x*f**(-slope[1]/g),y=self.y*f**(slope[0]/g))
            X = taylor(X,(self.x,0),(self.y,0),k-1)
        return X

    def _add_walls(self, diagram, k):
        P = self._path_product(diagram, k)
        for i in range(2,k-1):
            g = gcd(i-1, k-i-1)
            C = g*P.coefficient(self.x,i).coefficient(self.y,k-1-i)/(k-1-i)
            if C != 0:
                for j in range(2,len(diagram)):
                    wall = diagram[j]
                    slope = wall.slope
                    f = wall.function
                    if slope[0] == (i-1)/g and slope[1] == (k-1-i)/g:
                        replacement_wall = ScatteringWall( slope, f+C*self.x**(i-1)*self.y**(k-1-i) )
                        diagram[j] = replacement_wall 
                        break
                    elif slope[0] == 0 or slope[1]/slope[0] > RR(k-1-i)/(i-1):
                        new_wall = ScatteringWall( ((i-1)/g,(k-1-i)/g), 1+C*self.x**(i-1)*self.y**(k-1-i) )
                        diagram.insert(j, new_wall)
                        break
                    else:
                        j += 1
        return diagram

class BrokenLine(SageObject):
    """r
        A broken line
    """
    
    def __init__(self, initial_momentum):
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

#    #scattering methods
#    def scatter_monomial(self,momentum,alpha,beta,scatter_poly=-1):
#        #############
#        # Compute bad lands scattering polynomials here
#        #############
#
#        if scatter_poly == -1:
#            scatter_poly = 1 + self._x^(self._b*alpha)*self._y^(self._c*beta)
#        momentum_exp = vector(momentum.numerator().exponents()[0])-vector(momentum.denominator().exponents()[0])
#        
#        exp = abs(-beta*momentum_exp[0]+alpha*momentum_exp[1])
#        scattered_poly = momentum*scatter_poly^exp
#        #print scattered_poly
#        denom = scattered_poly.denominator()
#        numer = scattered_poly.numerator()
#        
#        scattered_monomials = []
#        for mon in numer.monomials():
#            scattered_monomials.append([numer.monomial_coefficient(mon),mon/denom,(alpha,beta)])
#        
#        return scattered_monomials
#        
#        
#    def monomial_degree(self,monomial):
#        return vector(monomial.numerator().exponents()[0])-vector(monomial.denominator().exponents()[0])
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
#            root = (sqrt((sqrt(self._c/self._b)^(n%2)*chebyshev_U(n,sqrt(self._b*self._c)/2))^2),sqrt((sqrt(self._b/self._c)^((n-1)%2)*chebyshev_U(n-1,sqrt(self._b*self._c)/2))^2))
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
#            root = (sqrt((sqrt(self._c/self._b)^((n-1)%2)*chebyshev_U(n-1,sqrt(self._b*self._c)/2))^2),sqrt((sqrt(self._b/self._c)^(n%2)*chebyshev_U(n,sqrt(self._b*self._c)/2))^2))
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
#    def broken_lines(self,init_momentum,Q):
#        (clockwise_walls,counterclockwise_walls) = self.find_walls(init_momentum)
#
#        clockwise_lines = [[[1,init_momentum,(0,0)]]]
#        for wall in clockwise_walls:
#            temp_lines = []
#            for line in clockwise_lines:
#                temp_lines.append(line)
#                last_momentum = line[-1]
#                for momentum_pair in self.scatter_monomial(last_momentum[1],wall[0],wall[1]):
#                    if momentum_pair[1].denominator().exponents()[0][1]>0 and momentum_pair[1] != last_momentum[1]:
#                        temp_lines.append(copy(line)+[momentum_pair])
#            clockwise_lines = temp_lines
#
#        counterclockwise_lines = [[[1,init_momentum,(0,0)]]]
#        for wall in counterclockwise_walls:
#            temp_lines = []
#            for line in counterclockwise_lines:
#                temp_lines.append(line)
#                last_momentum = line[-1]
#                for momentum_pair in self.scatter_monomial(last_momentum[1],wall[0],wall[1]):
#                    if momentum_pair[1].denominator().exponents()[0][0]>0 and momentum_pair[1] != last_momentum[1]:
#                        temp_lines.append(copy(line)+[momentum_pair])
#            counterclockwise_lines = temp_lines
#        
#        broken_lines = []
#        for line in clockwise_lines:
#            final_direction = self.monomial_degree(line[-1][1])
#            if final_direction[0] > 0:
#                broken_lines.append(line)
#            elif Q[0]*final_direction[1]-Q[1]*final_direction[0] < 0:
#                broken_lines.append(line)
#        for line in counterclockwise_lines:
#            final_direction = self.monomial_degree(line[-1][1])
#            if final_direction[1] > 0:
#                broken_lines.append(line)
#            elif Q[0]*final_direction[1]-Q[1]*final_direction[0] > 0:
#                broken_lines.append(line)
#        
#        return broken_lines
#
#
#    def GreedyElementBrokenLine(self,a1,a2,Q):
#        """
#        Output:
#            -Returns the greedy basis element
#        """
#
#        output = 0
#        for line in self.broken_lines(self._x^(-a1)*self._y^(-a2),Q):
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
#        walls = self.find_walls(self._x^(-a1)*self._y^(-a2))
#        length = ceil(sqrt(Q[0]^2+Q[1]^2))
#        for w in walls[0][:-1]+walls[1][:-1]:
#            wall_length = ceil(sqrt(w[0]^2+w[1]^2))
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
#        #creates tikz pictures of all broken lines with initial momentum x^{-a1}y^{-a2} ending at Q
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
#                TeXFile.write("^{")
#                if a1 > 0:
#                    TeXFile.write("-")
#                    TeXFile.write(str(a1)+"}")
#                else:
#                    TeXFile.write(str(-a1)+"}")
#        if a2 != 0:
#            TeXFile.write("y")
#            if a2 != -1:
#                TeXFile.write("^{")
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
#        for line in self.broken_lines(self._x^(-a1)*self._y^(-a2),Q):
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
#
#    def _thetad(self, scatter_diagram, k):
#        X = self._x
#        for (slope,f) in scatter_diagram:
#            g = gcd(slope)
#            X = X(x=self._x*f^(-slope[1]/g),y=self._y*f^(slope[0]/g))
#            X = taylor(X,(self._x,0),(self._y,0),k-1)
#        return X
#
#    def _scatter(self, scatter_diagram, k):
#        P = self._thetad(scatter_diagram, k)
#        scatter_diagram = list(scatter_diagram)
#        for i in range(2,k-1):
#            g = gcd(i-1, k-i-1)
#            C = g*P.coefficient(self._x,i).coefficient(self._y,k-1-i)/(k-1-i)
#            for j in range(2,len(scatter_diagram)):
#                (slope,f) = scatter_diagram[j]
#                if C == 0:
#                    break
#                if slope[0] == (i-1)/g and slope[1] == (k-1-i)/g:
#                    new_ray = (slope, f+C*self._x^(i-1)*self._y^(k-1-i))
#                    C = 0
#                    scatter_diagram[j] = new_ray 
#                elif slope[0] == 0 or slope[1]/slope[0] >= (k-1-i)/(i-1):
#                    new_ray = (((i-1)/g,(k-1-i)/g),1+C*self._x^(i-1)*self._y^(k-1-i))
#                    scatter_diagram.insert(j, new_ray)
#                    C = 0
#                else:
#                    j += 1
#        return tuple(scatter_diagram)
#
#    def _get_scattering_diagram(self,fx, fy, ktop):
#        scattering_diagram = (((-1,0),fx),((0,-1),fy),((1,0),fx),((0,1),fy))
#        for k in range(2,ktop+1):
#           scattering_diagram = self._scatter(scattering_diagram,k)
#        return scattering_diagram
#
