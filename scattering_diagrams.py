from sage.misc.cachefunc import cached_method                                                                                                                             

class ScatteringRing(SageObject):
    """
        A rank-2 scattering ring
        
        b: exponent of x
        c: exponent of y
        end_point: a point in the first quadrant

        EXAMPLES:

        sage: R=ScatteringRing(2,2); R
        The rank-2 scattering ring associated to the integers (2, 2).
        sage: R=ScatteringRing(2,2,end_point=(10,3));
        The rank-2 scattering ring associated to the integers (2, 2), with end
        point (10, 3).
    """
    def __init__(self, b, c, end_point=None):
        self._b = b
        self._c = c
        self.scatter_ring = LaurentPolynomialRing(QQ, 'x,y')
        self.scatter_field = self.scatter_ring.fraction_field()
        (self.x, self.y) = self.scatter_ring.gens()
        self._end_point=end_point
   
    def __repr__(self):
        output = "The rank-2 scattering ring associated to the integers " 
        output += str((self._b,self._c))
        if self._end_point is not None:
            output += ", with end point "
            output += str(self._end_point)
        output +="."
        return output

    @cached_method
    def scattering_diagram(self, max_degree):
        """
            Return all the walls of the scattering diagram associated to self
            with monomials of totat degree at most max_degree

            EXAMPLES:

            sage: R.scattering_diagram(16)
            A scattering diagram with 11 walls in it.
        """

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
    def broken_lines(self, init_momentum, end_point=None):
        """
        Return all the broken lines with initial momentum init_momentum and end
        point end_point.

        EXAMPLES:

        sage: R.broken_lines(R.x^-1*R.y^-1,(10,3))
        A collection of 3 broken lines with initial momentum x^-1*y^-1 and end
        point (10, 3).
        """
        if end_point is None:
            if self._end_point is not None:
                end_point = self._end_point
            else: 
                raise ValueError("You should specify the end point")

        degree = _monomial_degree(init_momentum)
        if all( x >= 0 for x in degree ):
            return BrokenLines(init_momentum, end_point, (BrokenLine(init_momentum),))
        elif degree[0] > 0:
            total_degree = degree[0]-degree[1]
            clockwise_diagram = ScatteringDiagram([ScatteringWall((-1,0),1+self.x**self._b)])
            counterclockwise_diagram = ScatteringDiagram([])
        elif degree[1] > 0:
            total_degree = -degree[0]+degree[1]
            clockwise_diagram = ScatteringDiagram([])
            counterclockwise_diagram = ScatteringDiagram([ScatteringWall((0,-1),1+self.x**self._c)])
        else:
            total_degree = -sum(degree)
            diagram = self.scattering_diagram(total_degree)
        
            clockwise_diagram = ScatteringDiagram([ W for W in diagram if _side(W.slope,degree) == "clockwise"  ])
            counterclockwise_diagram = ScatteringDiagram(reversed([ W for W in diagram if _side(W.slope,degree) == "counterclockwise"  ]))

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
        
        broken_lines = BrokenLines(init_momentum, end_point)
        for line in clockwise_lines:
            final_direction = _monomial_degree(line.monomial())
            if final_direction[0] > 0:
                broken_lines.add_line(line)
            elif end_point[0]*final_direction[1]-end_point[1]*final_direction[0] < 0:
                broken_lines.add_line(line)
        for line in counterclockwise_lines:
            final_direction = _monomial_degree(line.monomial())
            if final_direction[1] > 0:
                broken_lines.add_line(line)
            elif end_point[0]*final_direction[1]-end_point[1]*final_direction[0] > 0:
                broken_lines.add_line(line)
        
        return broken_lines
    
    def _scatter_at_wall(self, momentum, wall):
        momentum_exp = _monomial_degree(momentum)
        exp = abs(wall.slope[0]*momentum_exp[1]-wall.slope[1]*momentum_exp[0])
        scattered_function = momentum*(wall.function**exp)
        scattered_monomials = scattered_function.monomials()
        scattered_segments = []
        for mon in scattered_monomials:
            scattered_segments.append(BrokenLineSegment(scattered_function.monomial_coefficient(mon), mon, wall.slope))
        return tuple(scattered_segments)

    def show_broken_lines(self, initial_momentum, end_point=None, bounding_box=((-30,-30),(30,30))):
        if end_point is None:
            if self._end_point is not None:
                end_point = self._end_point
            else: 
                raise ValueError("You should specify the end point")

        tikz_commands = "\clip " +str(bounding_box[0]) +" rectangle " +str(bounding_box[1]) +";\n"
        tikz_commands += self.scattering_diagram(-sum(_monomial_degree(initial_momentum))).tikz(bounding_box=bounding_box)
        tikz_commands += self.broken_lines(initial_momentum, end_point).tikz(bounding_box=bounding_box)
        _show(tikz_commands)

    def save_broken_lines(self, initial_momentum, end_point=None, filename=None, bounding_box=((-30,-30),(30,30))):
        if end_point is None:
            if self._end_point is not None:
                end_point = self._end_point
            else: 
                raise ValueError("You should specify the end point")

        if filename is None:
            raise ValueError("A filename should be provided")
        tikz_commands = "\clip " +str(bounding_box[0]) +" rectangle " +str(bounding_box[1]) +";\n"
        tikz_commands += self.scattering_diagram(-sum(_monomial_degree(initial_momentum))).tikz(bounding_box=bounding_box)
        tikz_commands += self.broken_lines(initial_momentum, end_point).tikz(bounding_box=bounding_box)
        _save(tikz_commands, filename=filename)

    def greedy_element(self, initial_momentun, end_point=None):
        """
            Returns the greedy element computed with broken lines

            initial_momentum: a monomial whose degree is the $d$-vector of the
                greedy element

            end_point: any point in the first quadrant

            EXAMPLES:

            sage: R.greedy_element(R.x^-1*R.y^-1,(10,3))
            x*y^-1 + x^-1*y + x^-1*y^-1

        """
        if end_point is None:
            if self._end_point is not None:
                end_point = self._end_point
            else: 
                raise ValueError("You should specify the end point")

        output = 0
        for line in self.broken_lines(initial_momentun, end_point):
            output += line.coefficient()*line.monomial()
        return output

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

    def tikz(self, picture_radius=10, show_label=False):
        output = "\\draw[color=black,line width=1pt] (0,0) -- "
        endpoint = -N(picture_radius * vector(self.slope)/norm(vector(self.slope)))
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

    def __iter__(self):
        return iter(self.walls)
    
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

    def tikz(self, show_label=False, bounding_box=((-30,-30),(30,30))):
        picture_radius = max( [bounding_box[1][0]-bounding_box[0][0],bounding_box[1][1]-bounding_box[0][1]] )
        output = ""
        for W in self.walls:
            output += W.tikz(picture_radius=picture_radius, show_label=show_label)
        return output
    
    def show(self, bounding_box=((-30,-30),(30,30))):
        tikz_commands = "\clip " +str(bounding_box[0]) +" rectangle " +str(bounding_box[1]) +";\n"
        tikz_commands += self.tikz(bounding_box=bounding_box)
        _show(tikz_commands)

    def save(self, filename=None, bounding_box=((-30,-30),(30,30))):
        if filename is None:
            raise ValueError("A filename should be provided")
        tikz_commands = "\clip " +str(bounding_box[0]) +" rectangle " +str(bounding_box[1]) +";\n"
        tikz_commands += self.tikz(bounding_box=bounding_box)
        _save(tikz_commands, filename=filename)

class BrokenLines(SageObject):
    """r
        A collection of broken lines through a point
    """

    def __init__(self, initial_momentum, end_point, broken_lines=()):
        self.initial_momentum = initial_momentum
        self.end_point = end_point
        self.broken_lines = broken_lines

    def __iter__(self):
        return iter(self.broken_lines)
    
    def __repr__(self):
        output = "A collection of "
        output += str(len(self.broken_lines))
        output += " broken line" 
        if len(self.broken_lines)>1:
            output +="s"
        output += " with initial momentum "
        output += str(self.initial_momentum)
        output += " and end point "
        output += str(self.end_point)
        output += "."
        return output

    def add_line(self, broken_line):
        self.broken_lines += (broken_line,)

    def tikz(self, bounding_box=((-30,-30),(30,30))):
        picture_radius = max( [bounding_box[1][0]-bounding_box[0][0],bounding_box[1][1]-bounding_box[0][1]] )
        output = ""
        for bl in self.broken_lines:
            output += bl.tikz(self.end_point, picture_radius)
        output += "\\draw[color=black,fill=black] ("+str(self.end_point[0])+","+str(self.end_point[1])+") circle (3pt);\n"
        return output

    def show(self, bounding_box=((-30,-30),(30,30))):
        tikz_commands = "\clip " +str(bounding_box[0]) +" rectangle " +str(bounding_box[1]) +";\n"
        tikz_commands += self.tikz(bounding_box=bounding_box)
        _show(tikz_commands)

    def save(self, filename=None, bounding_box=((-30,-30),(30,30))):
        if filename is None:
            raise ValueError("A filename should be provided")
        tikz_commands = "\clip " +str(bounding_box[0]) +" rectangle " +str(bounding_box[1]) +";\n"
        tikz_commands += self.tikz(bounding_box=bounding_box)
        _save(tikz_commands, filename=filename)


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

    def __iter__(self):
        return iter(self.line_segments)
    
    def append_segment(self, segment):
        self.line_segments += (segment,)

    def coefficient(self):
        return prod([ x.coefficient for x in self.line_segments ])

    def monomial(self):
        return self.line_segments[-1].monomial

    def tikz(self, end_point, picture_radius):
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
                final_point = (RR(current_point[0]+sqrt(picture_radius)*current_direction[0]),RR(current_point[1]+sqrt(picture_radius)*current_direction[1]))
            output += "\\draw[color=red,line width=1pt] ("
            output += str(current_point[0])+","+str(current_point[1])+") -- ("
            output += str(final_point[0])+","+str(final_point[1])+");\n"
            current_point = final_point
            if intersection_slope != None:
                output += "\\draw[color=red,fill=red] ("+str(final_point[0])+","+str(final_point[1])+") circle (3pt);\n"
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

def _show(tikz_commands):
    from sage.misc.temporary_file import tmp_filename
    filename = tmp_filename(ext=".pdf")
    _save(tikz_commands, filename)
    os.system('%s %s 2>/dev/null 1>/dev/null &' % (sage.misc.viewer.pdf_viewer(), filename))

def _save(tikz_commands, filename=None):
    ext = os.path.splitext(filename)[1].lower()
    if ext not in ['.pdf', '.tex']:
        raise ValueError("The extension must be either pdf o tex")
    (cwd, cfn) = os.path.split(filename)
    basename = cfn.rsplit(".")[0]
    TeXFile=open(cwd+"/"+basename+".tex",'w')
    TeXFile.write("\\documentclass[tikz,border=10pt]{standalone}\n")
    TeXFile.write("\\begin{document}\n\n")
    TeXFile.write("\\begin{tikzpicture}\n")
    TeXFile.write(tikz_commands)
    TeXFile.write("\\end{tikzpicture}\n\n")
    TeXFile.write("\\end{document}")
    TeXFile.close()
    if ext == ".pdf":
        import subprocess
        if subprocess.call(['pdflatex', '-halt-on-error', basename+".tex"], cwd=cwd, stdout=subprocess.PIPE) != 0:
            raise RuntimeError("Unable to compile " +str(basename) +".tex. Try to run\n $cd " +str(cwd) +"; pdflatex " +str(basename) +".tex")
