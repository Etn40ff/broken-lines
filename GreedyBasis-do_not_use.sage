class GreedyBasis2(SageObject):
    """
    This class computes the non-commutative greedy elements generalizing Lee-Li-Zelevinsky Dyck path model.
    
    Attributes:
        ._b
        ._c
        ._qbase_ring
        ._torus
        ._qtorus
        ._qX
        ._qX_inv
        ._standard_gens
        ._standard_qgens
    Methods:
        .__init__
        .character
        .r
        .U
    """
    
    def __init__(self, b, c):
        self._b = b
        self._c = c
        
        R.<x,y> = ZZ[]
        scatter_ring = R.fraction_field()
        self._x = scatter_ring(x)
        self._y = scatter_ring(y)
        S.<v> = ZZ[]
        Q = S.fraction_field()

        T.<x1,x2> = LaurentPolynomialRing(QQ,2)
        self._torus = T
        self._standard_gens = [x1^(-1)*x2^self._c+x1^(-1),x1,x2,x1^self._b*x2^(-1)+x2^(-1)]


    def GreedyCoeff(self,a1,a2,p,q):
        p = Integer(p)
        q = Integer(q)
        if p == 0 and q == 0:
            return 1
        sum1 = 0
        for k in range(1,p+1):
            bin = 0
            if a2-self._c*q+k-1 >= k:
                bin = binomial(a2-self._c*q+k-1,k)
            sum1 += (-1)^(k-1)*self.GreedyCoeff(a1,a2,p-k,q)*bin
        sum2 = 0
        for l in range(1,q+1):
            bin = 0
            if a1-self._b*p+l-1 >= l:
                bin = binomial(a1-self._b*p+l-1,l)
            sum2 += (-1)^(l-1)*self.GreedyCoeff(a1,a2,p,q-l)*bin
        #print "sum1=",sum1,"sum2=",sum2
        return max(sum1,sum2)


    def GreedyElement(self,a1,a2):
        """
        Output:
            -Returns the greedy basis element
        """

        Dyck_path = self.get_Dyck_path(a1,a2)
        horiz_edges = Dyck_path[1]
        vert_edges = Dyck_path[2]

        x1 = self._torus._gens[0]
        x2 = self._torus._gens[1]
        output = 0
        for horiz in self.power_set(range(1,a1+1)):
            for vert in self.power_set(range(1,a2+1)):
                    if self.is_compatible(a1,a2,horiz,vert):
                        output += x1^(self._b*vert.__len__()-a1)*x2^(self._c*horiz.__len__()-a2)
        return output 


    #Useful list methods
    def intersect(self,list1,list2):
        """
        returns the list which contains all elements in common to list1 and list2
        """
        intersection = []
        if list1.__len__() < list2.__len__():
            for i in list1:
                if list2.count(i) != 0:
                    intersection.append(i)
        else:
            for i in list2:
                if list1.count(i) != 0:
                    intersection.append(i)
        return intersection


    def compliment(self, list, sublist):
        """
        Input:
            -A list and a sublist
        Output:
            -The elements of list which are not contained in sublist
        """
        output = []
        for item in list:
            if sublist.count(item) == 0:
                output.append(item)
        return output


    def power_set(self,list):
        """
        returns a list containing all possible sublist of the given list
        """
        sublists = [[]]
        for i in range(0,list.__len__()):
            for j in range(0,2**i):
                sublist = sublists.pop(0)
                sublists.append(copy(sublist))
                sublist.append(list[i])
                sublists.append(copy(sublist))
        return sublists


    #Dyck path methods
    def get_Dyck_path(self,a1,a2):
        """
        Input:
            -Non-negative integers a1 and a2
        Output:
            -The set of vertices in the corresponding maximal Dyck path
        """
        
        a1 = Integer(a1) #fixes a type checking issue
        a2 = Integer(a2)
        if a1 != 0:
            Dyck_path_slope = a2/a1
        else:
            Dyck_path_slope = 1000000000
        Dyck_vertices = []
        pt = [0,0]
        Dyck_vertices.append(copy(pt))
        vert_subs = []
        horiz_subs = []
        while pt[0] != a1 or pt[1] != a2:
            if pt[0] == 0 and pt[0] != a1:
                pt[0] += 1
                Dyck_vertices.append(copy(pt))
                horiz_subs.append(len(Dyck_vertices)-1)
            elif pt[0] == a1 or (pt[1]+1)/pt[0] <= Dyck_path_slope:
                pt[1] += 1
                Dyck_vertices.append(copy(pt))
                vert_subs.append(len(Dyck_vertices)-1)
            else:
                pt[0] += 1
                Dyck_vertices.append(copy(pt))
                horiz_subs.append(len(Dyck_vertices)-1)
        return [Dyck_vertices,horiz_subs,vert_subs]
        
        
    def horiz_indices(self,a1,a2,H):
        """
        Input:
            -subset H of the horizontal edges
        Output:
            -corresponding subset of the labels 1,...,a1
        """
        
        Dyck_path = self.get_Dyck_path(a1,a2)
        horiz_edges = Dyck_path[1]
        
        indices = []
        for h in H:
            indices.append(horiz_edges.index(h)+1)
        return indices
        
        
    def vert_indices(self,a1,a2,V):
        """
        Input:
            -subset V of the vertical edges
        Output:
            -corresponding subset of the labels 1,...,a2
        """
        
        Dyck_path = self.get_Dyck_path(a1,a2)
        vert_edges = Dyck_path[2]
        
        indices = []
        for v in V:
            indices.append(vert_edges.index(v)+1)
        return indices


    def local_horiz_shadows(self,a1,a2,H):
        """
        Input:
            -subset H of the horizontal edges labeled by 1,...,a1
        Output:
            -a dictionary whose labels are elements h of H, where the entry associated
             to h is the subset of vertical edges contained in the local shadow of h
        """

        Dyck_path = self.get_Dyck_path(a1,a2)
        horiz_edges = Dyck_path[1]
        vert_edges = Dyck_path[2]

        shadow_dict = dict()
        for h in H:
            local_shadow = []
            shadow_height = self._c
            position = horiz_edges[h-1]+1
            total_path_switch = false
            while local_shadow.__len__() < shadow_height and not total_path_switch:
                if position > a1+a2:
                    position = position-a1-a2
                if vert_edges.count(position) != 0:
                    local_shadow.append(vert_edges.index(position)+1)
                elif H.count(horiz_edges.index(position)+1) != 0:
                    shadow_height += self._c
                position += 1
                if position == horiz_edges[h-1]+1:
                    total_path_switch = true
            shadow_dict.setdefault(h,local_shadow)
        return shadow_dict


    def horiz_shadow(self,a1,a2,H):
        """
        Input:
            -subset H of the horizontal edges labeled by 1,...,a1
        Output:
            -the subset of vertical edges contained in the shadow of H
        """
        shadow = []
        local_horiz_shadows = self.local_horiz_shadows(a1,a2,H)
        for h in H:
            local_shadow = local_horiz_shadows[h]
            for v in local_shadow:
                if shadow.count(v) == 0:
                    shadow.append(v)
            shadow.sort()
        return shadow
        
        
    def remote_horiz_shadow(self,a1,a2,H):
        """
        Input:
            -subset H of the horizontal edges labeled by 1,...,a1
        Output:
            -the subset of vertical edges contained in the remote shadow of H
             i.e. those vertical edges which potentially could be compatible with H
        """
        horiz_shadow = self.horiz_shadow(a1,a2,H)
        remote_shadow = []
        for v in horiz_shadow:
           if self.is_compatible(a1,a2,H,[v]):
               remote_shadow.append(v)
        return remote_shadow


    def vert_shadow(self,a1,a2,V):
        """
        Input:
            -subset V of the vertical edges labeled by 1,...,a2
        Output:
            -the subset of horizontal edges contained in the shadow of V
        """

        Dyck_path = self.get_Dyck_path(a1,a2)
        horiz_edges = Dyck_path[1]
        vert_edges = Dyck_path[2]

        shadow = []
        if a1 == 0:
            return shadow
        for v in V:
            position = v - 1
            while vert_edges.count(position) != 0 and position > 0:
                position -= 1
            end = horiz_edges.index(position) + 1
            if end > self._b:
                start = end - self._b
            else:
                start = 0
            switch = true
            while shadow.count(horiz_edges[start]) != 0 and start != 0:
                if switch:
                    for i in range(start,end):
                        if shadow.count(horiz_edges[i]) == 0:
                            shadow.append(horiz_edges[i])
                switch = false
                end = start
                if end > self._b:
                    start = end - self._b
                else:
                    start = 0 
            for i in range(start,end):
                if shadow.count(horiz_edges[i]) == 0:
                    shadow.append(horiz_edges[i])
        shadow.sort()
        return shadow


    def is_compatible(self,a1,a2,H,V):
        """
        determines whether the subset H of the horizontal edges and the subset V of the vertical edges are compatible
        """

        Dyck_path = self.get_Dyck_path(a1,a2)
        horiz_edges = Dyck_path[1]
        vert_edges = Dyck_path[2]

        for h in H:
            for v in V:
                compatible = false
                p1 = horiz_edges[h-1]
                p2 = vert_edges[v-1]
                if p1 < p2:
                    for e in range(p1,p2+1):
                        #print "e=", e, range(p1,e+1), range(e,p2+1)
                        he1 = self.horiz_indices(a1,a2,self.intersect(range(p1,e+1),horiz_edges))
                        he2 = self.vert_indices(a1,a2,self.intersect(range(p1,e+1),vert_edges))
                        ev1 = self.horiz_indices(a1,a2,self.intersect(range(e,p2+1),horiz_edges))
                        ev2 = self.vert_indices(a1,a2,self.intersect(range(e,p2+1),vert_edges))
                        #print "he1=", he1
                        #print "he2=", he2
                        #print "ev1=", ev1
                        #print "ev2=", ev2
                        if he2.__len__() == self._c*self.intersect(he1,H).__len__() and e != p2:
                            #print "horizontal check"
                            compatible = true
                        if ev1.__len__() == self._b*self.intersect(ev2,V).__len__() and e != p1:
                            #print self.intersect(ev2,V)
                            #print "vertical check"
                            compatible = true
                        #print "Compatible?", compatible
                else:
                    for e in range(p1,a1+a2+p2+1):
                        if e <= a1+a2:
                            R1 = range(p1,e+1)
                            R2 = range(e,a1+a2+1)
                            R2.extend(range(1,p2+1))
                        else:
                            R1 = range(p1,a1+a2+1)
                            R1.extend(range(1,e-a1-a2+1))
                            R2 = range(e-a1-a2,p2+1)
                        #print "e=", e
                        #print "R1=", R1
                        #print "R2=", R2
                        he1 = self.horiz_indices(a1,a2,self.intersect(R1,horiz_edges))
                        he2 = self.vert_indices(a1,a2,self.intersect(R1,vert_edges))
                        ev1 = self.horiz_indices(a1,a2,self.intersect(R2,horiz_edges))
                        ev2 = self.vert_indices(a1,a2,self.intersect(R2,vert_edges))
                        #print "he1=", he1
                        #print "he2=", he2
                        #print "ev1=", ev1
                        #print "ev2=", ev2
                        if he2.__len__() == self._c*self.intersect(he1,H).__len__() and e != p2:
                            compatible = true
                        if ev1.__len__() == self._b*self.intersect(ev2,V).__len__() and e != p1:
                            compatible = true
                        #print "Compatible?", compatible
                if not compatible:
                    return false
        return true
        
        
    def get_restricted_compatible_pairs(self,a1,a2,set_horiz_card=-1,set_vert_card=-1):
        """
        Output:
            -returns compatible pairs (H,V) where V is in the shadow of H
        """

        Dyck_path = self.get_Dyck_path(a1,a2)
        horiz_edges = Dyck_path[1]
        vert_edges = Dyck_path[2]
        
        compatible_pairs = []
        for H in self.power_set(range(1,a1+1)):
            if set_horiz_card == -1 or H.__len__() == set_horiz_card:
                remote_shadow = self.remote_horiz_shadow(a1,a2,H)
                for V in self.power_set(remote_shadow):
                    if set_vert_card == -1 or V.__len__() == set_vert_card:
                        if self.is_compatible(a1,a2,H,V):
                            compatible_pairs.append([H,V])
        return compatible_pairs


    def GreedyElementRecursive(self,a1,a2):
        x1 = self._torus._gens[0]
        x2 = self._torus._gens[1]
        if a1 < 0:
            if a2 < 0:
                return x1^(-a1)*x2^(-a2)
            else:
                return x1^(-a1)*self._standard_gens[3]^a2
        elif a2 < 0:
            return self._standard_gens[0]^a1*x2^(-a2)
        output = 0
        for p in range(0,a2+1):
            for q in range(0,a1+1):
                output += self.GreedyCoeff(a1,a2,p,q)*x1^(self._b*p)*x2^(self._c*q)
        return x1^(-a1)*x2^(-a2)*output


    def GreedyTest(self):
        for a1 in range(1,6):
            for a2 in range(1,6):
                print "a1=",a1,"a2=",a2,self.GreedyElement(a1,a2) == self.GreedyElementRecursive(a1,a2)
        

    def draw_diagonal(self,a1,a2):
        #produces the tikz string for drawing the main diagonal of the rectangle R_self._n
        output = "  \\draw[step=0.25cm,color=gray] (0,0) grid ("+str(float(a1/4))+","+str(float(a2/4))+");\n"
        output += "  \\draw[color=gray] (0,0) -- ("+str(float(a1/4))+","+str(float(a2/4))+");\n"
        return output


    def draw_Dyck_path_vertices(self,a1,a2):
        #produces the tikz string for drawing the vertices of the maximal Dyck path
        output = ""
        Dyck_vertices = self.get_Dyck_path(a1,a2)[0]
        for vert in Dyck_vertices:
            output += "  \\draw[fill=black] ("+str(float(vert[0]/4))+","+str(float(vert[1]/4))+") circle (1.1pt);\n"
        return output
        
        
    def draw_compatible_pair(self,a1,a2,cp):
        
        H = cp[0]
        V = cp[1]
        
        Dyck_path = self.get_Dyck_path(a1,a2)
        Dyck_vertices = Dyck_path[0]
        horiz_edges = Dyck_path[1]
        vert_edges = Dyck_path[2]
        
        output = self.draw_diagonal(a1,a2)
        
        for h in H:
            vertex1 = Dyck_vertices[horiz_edges[h-1]-1]
            vertex2 = Dyck_vertices[horiz_edges[h-1]]
            output += "  \\draw[color=red,line width=1.5pt] ("
            output += str(float(vertex1[0]/4))+","+str(float(vertex1[1]/4))+") -- ("
            output += str(float(vertex2[0]/4))+","+str(float(vertex2[1]/4))+");\n"
        
        for v in V:
            vertex1 = Dyck_vertices[vert_edges[v-1]-1]
            vertex2 = Dyck_vertices[vert_edges[v-1]]
            output += "  \\draw[color=red,line width=1.5pt] ("
            output += str(float(vertex1[0]/4))+","+str(float(vertex1[1]/4))+") -- ("
            output += str(float(vertex2[0]/4))+","+str(float(vertex2[1]/4))+");\n"
        
        for v in self.compliment(range(1,a2+1),self.horiz_shadow(a1,a2,H)):
            vertex1 = Dyck_vertices[vert_edges[v-1]-1]
            vertex2 = Dyck_vertices[vert_edges[v-1]]
            output += "  \\draw[color=green,line width=1.5pt] ("
            output += str(float(vertex1[0]/4))+","+str(float(vertex1[1]/4))+") -- ("
            output += str(float(vertex2[0]/4))+","+str(float(vertex2[1]/4))+");\n"
            
        output += self.draw_Dyck_path_vertices(a1,a2)
            
        return output
        
    
    def draw_compatible_pairs(self, a1, a2, max_counter=1,set_horiz_card=-1,set_vert_card=-1):
        #creates tikz pictures of all compatible pairs on the maximal Dyck path in the rectangle (a1,a2)
        #gray edges may be freely included

        working_dir="/tmp/"
        #filename = "comp_pairs"+str(self._b)+str(self._c)+"-("+str(a1)+","+str(a2)+")"
        #if set_horiz_card != -1:
        #    filename += "-|H|="+str(set_horiz_card)
        #if set_vert_card != -1:
        #    filename += "-|V|="+str(set_vert_card)
        filename = "sage_output"
        filename += ".tex"
        print filename
        TeXFile=open(working_dir+filename,'w')
        TeXFile.write("\\documentclass{article}\n")
        TeXFile.write("\\usepackage{amsmath, amssymb, latexsym, tikz}\n")
        TeXFile.write("\\usepackage{pgflibraryarrows,pgflibrarysnakes}\n\n")
        TeXFile.write("\\begin{document}\n\n")
        TeXFile.write("Let $b="+str(self._b)+"$ and $c="+str(self._c)+"$.  ")
        TeXFile.write("We consider compatible pairs in the Dyck path $D_{"+str(a1)+","+str(a2)+"}$ given by:\\\\\n\n")
        TeXFile.write(" \\begin{tikzpicture}\n")
        TeXFile.write(self.draw_diagonal(a1,a2)+self.draw_Dyck_path_vertices(a1,a2))
        TeXFile.write(" \\end{tikzpicture}\\\\\n\n")
        TeXFile.write("where a green edge may be freely included or excluded.  ")
        TeXFile.write("These compute the greedy element: $x["+str(a1)+","+str(a2)+"]=$\\\\\n\n")
        
        counter = 0
        for compatible_pair in self.get_restricted_compatible_pairs(a1,a2,set_horiz_card,set_vert_card):
            TeXFile.write(" \\begin{tikzpicture}\n")
            TeXFile.write(self.draw_compatible_pair(a1,a2,compatible_pair))
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
