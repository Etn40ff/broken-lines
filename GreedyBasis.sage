class GreedyBasis(SageObject):
    """r
        Class to manipulate Greedy Elements and Boken Lines
    """
    def __init__(self, b, c):
        self._b = b
        self._c = c
        
        var('x,y')
        self._scatter_ring=ZZ[x,y].fraction_field()
        self._x = self._scatter_ring(x)
        self._y = self._scatter_ring(y)
    
    def get_Dyck_path(self,end):
        """
            Maximal Dyck path from (0,0) to end
        """
        
        if any( x < 0 for x in end):
            raise ValueError("end point must be in the first quadrant")
        vertices=[(0,0)]
        horizontals=[]
        verticals=[]

        if end[0] == 0:
            for i in range(1,end[1]+1):
                p = (0,i)
                vertices.append(p)
                verticals.append(i)
        
        elif end[1] == 0:
            for i in range(1,end[0]+1):
                p = (i,0)
                vertices.append(p)
                horizontals.append(i)

        else:
            slope = end[1]/end[0]
            p = (1,0)
            vertices.append(p)
            horizontals.append(1)
            counter = 2
            while p != end:
                if (p[1]+1)/p[0] <= slope:
                    p=(p[0],p[1]+1)
                    verticals.append(counter)
                else:
                    p=(Integer(p[0]+1),Integer(p[1]))
                    horizontals.append(counter)
                vertices.append((p[0],p[1]))
                counter += 1
        return (tuple(vertices),tuple(horizontals),tuple(verticals))

    def _draw_grid(self,end):
        output = "  \\draw[step=0.25cm,color=gray] (0,0) grid ("+str(float(end[0]/4))+","+str(float(end[1]/4))+");\n"
        output += "  \\draw[color=gray] (0,0) -- ("+str(RR(end[0])/4)+","+str(RR(end[1])/4)+");\n"
        return output


    def _draw_Dyck_path_vertices(self,end):
        output = ""
        Dyck_vertices = self.get_Dyck_path(end)[0]
        for vert in Dyck_vertices:
            output += "  \\draw[fill=black] ("+str(RR(vert[0])/4)+","+str(RR(vert[1])/4)+") circle (1.1pt);\n"
        return output
    
    def _draw_compatible_pair(self,end,cp):
        H = cp[0]
        V = cp[1]
        dyck_path = self.get_Dyck_path(end)
        vertices = dyck_path[0]
        
        output = self._draw_grid(end)
        
        for h in H:
            vertex1 = vertices[h-1]
            vertex2 = vertices[h]
            output += "  \\draw[color=red,line width=1.5pt] ("
            output += str(RR(vertex1[0])/4)+","+str(RR(vertex1[1])/4)+") -- ("
            output += str(RR(vertex2[0])/4)+","+str(RR(vertex2[1])/4)+");\n"
        
        for v in V:
            vertex1 = vertices[v-1]
            vertex2 = vertices[v]
            output += "  \\draw[color=red,line width=1.5pt] ("
            output += str(RR(vertex1[0])/4)+","+str(RR(vertex1[1])/4)+") -- ("
            output += str(RR(vertex2[0])/4)+","+str(RR(vertex2[1])/4)+");\n"
            
        output += self._draw_Dyck_path_vertices(end)

        return output
    
    def latex_standalone_tikz(self, tikz_commands, working_dir="/tmp/", filename="sage_output"):
        import subprocess
        TeXFile=open(working_dir+filename+".tex",'w')
        TeXFile.write("\\documentclass[tikz,border=10pt]{standalone}\n")
        #TeXFile.write("\\usepackage{amsmath, amssymb, latexsym, tikz}\n")
        #TeXFile.write("\\usepackage{pgflibraryarrows,pgflibrarysnakes}\n\n")
        TeXFile.write("\\begin{document}\n\n")
        TeXFile.write("\\begin{tikzpicture}\n")
        TeXFile.write(tikz_commands)
        TeXFile.write("\\end{tikzpicture}\n\n")
        TeXFile.write("\\end{document}")
        TeXFile.close()
        subprocess.call(['pdflatex', '-halt-on-error', filename+".tex"], cwd=working_dir, stdout=subprocess.PIPE)

    def compatible_pairs(self, end):
        output=[]
        dyck_path = self.get_Dyck_path(end)
        pairs = [ ( [], [], list(dyck_path[1]), list(dyck_path[2]) ) ]
        while pairs:
            pair = pairs.pop()
            if pair[2]:
                x = pair[2].pop()
                pairs += [tuple(map(copy,pair))]
                if self.is_compatible_pair(( pair[0]+[x], pair[1], pair[2], pair[3] ), dyck_path):
                    pairs += [ ( pair[0]+[x], pair[1], pair[2], pair[3] ) ]
            elif pair[3]:
                x = pair[3].pop()
                pairs += [tuple(map(copy,pair))]
                if self.is_compatible_pair(( pair[0], pair[1]+[x], pair[2], pair[3] ), dyck_path):
                    pairs += [ ( pair[0], pair[1]+[x], pair[2], pair[3] ) ]
            else: 
                output += [tuple(map(tuple,pair[0:2]))]
        output.sort(key=lambda x: (len(x[0])+len(x[1]),len(x[0]),len(x[1])))
        return output

    def is_compatible_pair(self, pair, dyck_path):
        if len(pair[0])*len(pair[1]) == 0:
            return True
        for u in pair[0]:
            E = dyck_path[0][u-1]
            for v in pair[1]:
                F = dyck_path[0][v]
                subpath = self._subpath(E,F,dyck_path)
                found = False
                for A in [ p for p in subpath[0] if p not in [E,F] ]:
                    af = self._subpath(A,F,dyck_path)
                    if len(af[1]) == self._b*len([x for x in af[2] if x in pair[1]]):
                        found = True
                        break
                    ea = self._subpath(E,A,dyck_path)
                    if len(ea[2]) == self._c*len([x for x in ea[1] if x in pair[0]]):
                        found = True
                        break
                if not found:
                    return False
        return True
                
    def _subpath(self, p, q, dyck_path):
        p_index = dyck_path[0].index(p)
        q_index = dyck_path[0].index(q)
        if p_index == q_index:
            return dyck_path
        elif p_index < q_index:
            vertices = ( x for x in dyck_path[0] if self._follows(x,p) and
                    self._follows(q,x) )
            horizontals = (x for x in dyck_path[1] if x > p_index and 
                    x <= q_index )
            verticals = (x for x in dyck_path[2] if x > p_index and 
                    x <= q_index )
        else:
            vertices = ( x for x in dyck_path[0] if self._follows(x,p) or
                    self._follows(q,x) )
            horizontals = (x for x in dyck_path[1] if x > p_index or 
                    x <= q_index )
            verticals = (x for x in dyck_path[2] if x > p_index or 
                    x <= q_index )
        return map(tuple,(vertices, horizontals, verticals))

    def _follows(self,p,q):
        if p[0] >= q[0] and p[1] >= q[1]:
            return True
        return False

    def _appendable_edge(self, pair, dyck_path, t, edge):
        if t == "v":
            if pair[0] == []:
                return True
            else:
                F = dyck_path[0][edge]
                for i in pair[0]:
                    E = dyck_path[0][i-1]


        if t == "h" and pair[1] == []:
            return True
        
        return False

    def draw_compatible_pairs(self, end, max_counter=1, working_dir="/tmp/", filename="sage_output"):
        filename += ".tex"
        TeXFile=open(working_dir+filename,'w')
        TeXFile.write("\\documentclass{article}\n")
        TeXFile.write("\\usepackage{amsmath, amssymb, latexsym, tikz}\n")
        TeXFile.write("\\usepackage{pgflibraryarrows,pgflibrarysnakes}\n\n")
        TeXFile.write("\\begin{document}\n\n")
        TeXFile.write("Let $b="+str(self._b)+"$ and $c="+str(self._c)+"$.  ")
        TeXFile.write("We consider compatible pairs in the Dyck path $D_{"+str(end[0])+","+str(end[1])+"}$ given by:\\\\\n\n")
        TeXFile.write(" \\begin{tikzpicture}\n")
        TeXFile.write(self._draw_grid(end)+self._draw_Dyck_path_vertices(end))
        TeXFile.write(" \\end{tikzpicture}\\\\\n\n")
        TeXFile.write("These compute the greedy element: $x["+str(end[0])+","+str(end[1])+"]$.\n")
        TeXFile.write("(They are sorted lexicographically using (\#edges,\#horizontals,\#verticals) as key.)\\\\\n\n")
        
        counter = 0
        for compatible_pair in self.compatible_pairs(end):
            TeXFile.write(" \\begin{tikzpicture}\n")
            TeXFile.write(self._draw_compatible_pair(end,compatible_pair))
            TeXFile.write(" \\end{tikzpicture}\quad\n")
            counter += 1
            if counter == max_counter:
                counter = 0
                TeXFile.write("\n\\vspace{.3in}\n\n")
        
        TeXFile.write("\\end{document}")
        TeXFile.close()
        
        import subprocess
        subprocess.call(['pdflatex', '-halt-on-error', filename], cwd=working_dir, stdout=subprocess.PIPE)
        
    #scattering methods
    def scatter_monomial(self,momentum,alpha,beta,scatter_poly=-1):
        #############
        # Compute bad lands scattering polynomials here
        #############

        if scatter_poly == -1:
            scatter_poly = 1 + self._x^(self._b*alpha)*self._y^(self._c*beta)
        momentum_exp = vector(momentum.numerator().exponents()[0])-vector(momentum.denominator().exponents()[0])
        
        exp = abs(-beta*momentum_exp[0]+alpha*momentum_exp[1])
        scattered_poly = momentum*scatter_poly^exp
        #print scattered_poly
        denom = scattered_poly.denominator()
        numer = scattered_poly.numerator()
        
        scattered_monomials = []
        for mon in numer.monomials():
            scattered_monomials.append([numer.monomial_coefficient(mon),mon/denom,(alpha,beta)])
        
        return scattered_monomials
        
        
    def monomial_degree(self,monomial):
        return vector(monomial.numerator().exponents()[0])-vector(monomial.denominator().exponents()[0])


    def find_walls(self,init_momentum):
        counterclockwise_root_walls = [(0,1)]
        clockwise_root_walls = [(1,0)]
        counterclockwise_split_root_walls = []
        clockwise_split_root_walls = []
        
        momentum_exp = vector(init_momentum.numerator().exponents()[0])-vector(init_momentum.denominator().exponents()[0])
        check1 = true
        n=0
        while check1:
            root = (sqrt((sqrt(self._c/self._b)^(n%2)*chebyshev_U(n,sqrt(self._b*self._c)/2))^2),sqrt((sqrt(self._b/self._c)^((n-1)%2)*chebyshev_U(n-1,sqrt(self._b*self._c)/2))^2))
            check1 = false
            if momentum_exp[0]*root[1]-momentum_exp[1]*root[0]>0:
                if momentum_exp[0]+self._b*root[0]<0:
                    counterclockwise_root_walls = [root]+counterclockwise_root_walls
                    check1 = true
            elif momentum_exp[0]*root[1]-momentum_exp[1]*root[0]<0:
                if momentum_exp[1]+self._c*root[1]<0:
                    counterclockwise_split_root_walls.append(root)
                    check1 = true
            n += 1
        check2 = true
        n=0
        while check2:
            root = (sqrt((sqrt(self._c/self._b)^((n-1)%2)*chebyshev_U(n-1,sqrt(self._b*self._c)/2))^2),sqrt((sqrt(self._b/self._c)^(n%2)*chebyshev_U(n,sqrt(self._b*self._c)/2))^2))
            check2 = false
            if momentum_exp[0]*root[1]-momentum_exp[1]*root[0]<0:
                if momentum_exp[1]+self._c*root[1]<0:
                    clockwise_root_walls = [root]+clockwise_root_walls
                    check2 = true
            elif momentum_exp[0]*root[1]-momentum_exp[1]*root[0]>0:
                if momentum_exp[0]+self._b*root[0]<0:
                    clockwise_split_root_walls.append(root)
                    check2 = true
            n += 1

        ###########
        # Imaginary walls needed
        ###########
        left_imaginary_walls = []
        right_imaginary_walls = []
        counterclockwise_walls = clockwise_split_root_walls+left_imaginary_walls+counterclockwise_root_walls
        clockwise_walls = counterclockwise_split_root_walls+right_imaginary_walls+clockwise_root_walls

        #print "clockwise=",clockwise_walls
        #print "counterclockwise_walls=",counterclockwise_walls

        return [clockwise_walls,counterclockwise_walls]
        
        
    def broken_lines(self,init_momentum,Q):
        (clockwise_walls,counterclockwise_walls) = self.find_walls(init_momentum)

        clockwise_lines = [[[1,init_momentum,(0,0)]]]
        for wall in clockwise_walls:
            temp_lines = []
            for line in clockwise_lines:
                temp_lines.append(line)
                last_momentum = line[-1]
                for momentum_pair in self.scatter_monomial(last_momentum[1],wall[0],wall[1]):
                    if momentum_pair[1].denominator().exponents()[0][1]>0 and momentum_pair[1] != last_momentum[1]:
                        temp_lines.append(copy(line)+[momentum_pair])
            clockwise_lines = temp_lines

        counterclockwise_lines = [[[1,init_momentum,(0,0)]]]
        for wall in counterclockwise_walls:
            temp_lines = []
            for line in counterclockwise_lines:
                temp_lines.append(line)
                last_momentum = line[-1]
                for momentum_pair in self.scatter_monomial(last_momentum[1],wall[0],wall[1]):
                    if momentum_pair[1].denominator().exponents()[0][0]>0 and momentum_pair[1] != last_momentum[1]:
                        temp_lines.append(copy(line)+[momentum_pair])
            counterclockwise_lines = temp_lines
        
        broken_lines = []
        for line in clockwise_lines:
            final_direction = self.monomial_degree(line[-1][1])
            if final_direction[0] > 0:
                broken_lines.append(line)
            elif Q[0]*final_direction[1]-Q[1]*final_direction[0] < 0:
                broken_lines.append(line)
        for line in counterclockwise_lines:
            final_direction = self.monomial_degree(line[-1][1])
            if final_direction[1] > 0:
                broken_lines.append(line)
            elif Q[0]*final_direction[1]-Q[1]*final_direction[0] > 0:
                broken_lines.append(line)
        
        return broken_lines


    def GreedyElementBrokenLine(self,a1,a2,Q):
        """
        Output:
            -Returns the greedy basis element
        """

        output = 0
        for line in self.broken_lines(self._x^(-a1)*self._y^(-a2),Q):
            coeff = prod([i[0] for i in line])
            output += coeff*line[-1][1]
        return output 
        
    
    def draw_broken_line(self,line,Q,rgb=(1,0,0)):
        color = "{rgb:red,"+str(rgb[0])+";green,"+str(rgb[1])+";blue,"+str(rgb[2])+"}"
        output = ""
        current_point = Q
        for i in range(1,line.__len__()):
            current_direction  = self.monomial_degree(line[-i][1])
            intersection_slope = line[-i][2]
            #print "line=",line[-i]
            #print "direction=",current_direction
            #print "intersection=",intersection_slope
            var('s','t')
            solution_dict = solve([current_direction[1]*(s-current_point[0])==current_direction[0]*(t-current_point[1]),intersection_slope[1]*s==intersection_slope[0]*t],(s,t),solution_dict=True)[0]
            final_point = (solution_dict[s],solution_dict[t])
            output += "  \\draw[color="+color+",line width=1pt] ("
            output += str(current_point[0])+","+str(current_point[1])+") -- ("
            output += str(final_point[0])+","+str(final_point[1])+");\n"
            output += "  \\draw[color="+color+",fill="+color+"] ("+str(final_point[0])+","+str(final_point[1])+") circle (3pt);\n"
            current_point = final_point
        current_direction  = self.monomial_degree(line[0][1])
        final_point = (current_point[0]+max(Q)*current_direction[0],current_point[1]+max(Q)*current_direction[1])
        output += "  \\draw[color="+color+",line width=1pt] ("
        output += str(current_point[0])+","+str(current_point[1])+") -- ("
        output += str(final_point[0])+","+str(final_point[1])+");\n"
        return output


    def draw_walls(self,a1,a2,Q):
        output = ""
        walls = self.find_walls(self._x^(-a1)*self._y^(-a2))
        length = ceil(sqrt(Q[0]^2+Q[1]^2))
        for w in walls[0][:-1]+walls[1][:-1]:
            wall_length = ceil(sqrt(w[0]^2+w[1]^2))
            output += "  \\draw[color=black,line width=1pt] (0,0) -- ("
            output += str(-ceil(self._b*length*w[0]/wall_length))+","+str(-ceil(self._c*length*w[1]/wall_length))+");\n"
        for w in [(0,1),(1,0)]:
            output += "  \\draw[color=black,line width=1pt] (0,0) -- ("
            output += str(self._b*length*w[0])+","+str(self._c*length*w[1])+");\n"
        #print output
        return output
        
        
    def draw_broken_lines(self,a1,a2,Q):
        #creates tikz pictures of all broken lines with initial momentum x^{-a1}y^{-a2} ending at Q
        working_dir="/tmp/"
        #filename = "broken_lines"+str(self._b)+str(self._c)+"-("+str(a1)+","+str(a2)+")"#-("+str(Q[0])+","+str(Q[1])+")"
        filename = "sage_output"
        filename += ".tex"
        print filename
        
        TeXFile=open(working_dir+filename,'w')
        TeXFile.write("\\documentclass{article}\n")
        TeXFile.write("\\usepackage{amsmath, amssymb, latexsym, tikz}\n")
        TeXFile.write("\\usepackage{pgflibraryarrows,pgflibrarysnakes}\n\n")
        TeXFile.write("\\begin{document}\n\n")
        TeXFile.write("Let $b="+str(self._b)+"$ and $c="+str(self._c)+"$.  ")
        TeXFile.write("We consider broken lines with initial momentum $")
        if a1 != 0:
            TeXFile.write("x")
            if a1 != -1:
                TeXFile.write("^{")
                if a1 > 0:
                    TeXFile.write("-")
                    TeXFile.write(str(a1)+"}")
                else:
                    TeXFile.write(str(-a1)+"}")
        if a2 != 0:
            TeXFile.write("y")
            if a2 != -1:
                TeXFile.write("^{")
                if a2 > 0:
                    TeXFile.write("-")
                    TeXFile.write(str(a2)+"}")
                else:
                    TeXFile.write(str(-a2)+"}")
        
        TeXFile.write("$ which end at $Q=("+str(Q[0])+","+str(Q[1])+")$.  These are given by:\\\\\n\n")
        TeXFile.write("\\resizebox{4.5in}{4.5in}{\n")
        TeXFile.write(" \\begin{tikzpicture}\n")
        TeXFile.write(self.draw_walls(a1,a2,Q))
        
        for line in self.broken_lines(self._x^(-a1)*self._y^(-a2),Q):
            TeXFile.write(self.draw_broken_line(line,Q))


        TeXFile.write("  \\draw[fill=black] ("+str(Q[0])+","+str(Q[1])+") circle (3pt);\n")
        
        TeXFile.write(" \\end{tikzpicture}}\\\\\n\n")
            
        TeXFile.write("\\end{document}")
        TeXFile.close()
        
        import subprocess
        subprocess.call(['pdflatex', '-halt-on-error', filename], cwd=working_dir, stdout=subprocess.PIPE)
