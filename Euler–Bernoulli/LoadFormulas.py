'''
BSD 3-Clause License

Copyright (c) 2019, open-struct-engineer developers and Donald N. Bockoven III
All rights reserved.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:

* Redistributions of source code must retain the above copyright notice, this
  list of conditions and the following disclaimer.

* Redistributions in binary form must reproduce the above copyright notice,
  this list of conditions and the following disclaimer in the documentation
  and/or other materials provided with the distribution.

* Neither the name of the copyright holder nor the names of its
  contributors may be used to endorse or promote products derived from
  this software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
'''

from __future__ import division
from numpy import sign
from numpy import zeros
import numpy as np
import math

def PieceFunctionString(piece_set):
        '''
        # Returns the general piecwise function in the form of a string
        # INPUT: List
        # List makup:
        # list1 is the polynomial coeficients of order [c0,c1x,c2x^2,...,cnx^n]
        # where the list values will only by the cn's*
        # list 2 will be the range over which the function piece applies
        # 0 <= a would be [0,a] **note it will be assumed the the eqality is <= not <
        # rerturned lists will be [[[list11],[list21]],....,[[list1n],[list2n]]
        # where n is the total number of functions to capture the range from
        # 0 to the full span, L of the beam
        '''

        output = ''

        for func in piece_set:
            i=0

            if all(c == 0 for c in func[0]):
                line = '0'
            else:
                line = ''
                for c in func[0]:
                    if c == 0:
                        pass
                    elif i == 0:
                        line = line + '{0:0.4f}'.format(c)
                    elif i == 1:
                        if line == '':
                            line = line + '{0:0.4f}*x'.format(c)
                        elif c < 0:
                            line = line + '-{0:0.4f}*x'.format(abs(c))
                        else:
                            line = line + '+{0:0.4f}*x'.format(c)
                    else:
                        if line == '':
                            line = line + '{0:0.4f}*x^{1}'.format(c,i)
                        elif c < 0:
                            line = line + '-{0:0.4f}*x^{1}'.format(abs(c),i)
                        else:
                            line = line + '+{0:0.4f}*x^{1}'.format(c,i)
                    i+=1

            output = output + '{0:0.4f} < x <= {1:0.4f}:\n'.format(func[1][0],func[1][1]) + line + '\n'
        return output

def PieceFunctionStringHTMLTable(piece_set,heading_str):
    '''
    # Returns the general piecwise function in the form of a string
    # INPUT: List
    # List makup:
    # list1 is the polynomial coeficients of order [c0,c1x,c2x^2,...,cnx^n]
    # where the list values will only by the cn's*
    # list 2 will be the range over which the function piece applies
    # 0 <= a would be [0,a] **note it will be assumed the the eqality is <= not <
    # rerturned lists will be [[[list11],[list21]],....,[[list1n],[list2n]]
    # where n is the total number of functions to capture the range from
    # 0 to the full span, L of the beam
    '''

    output = '<table>\n<tr>\n<th>{0}</th>\n</tr>\n'.format(heading_str)

    for func in piece_set:
        i=0

        if all(c == 0 for c in func[0]):
            line = '0'
        else:
            line = ''
            for c in func[0]:
                if c == 0:
                    pass
                elif i == 0:
                    line = line + '{0:0.4f}'.format(c)
                elif i == 1:
                    if line == '':
                        line = line + '{0:0.4f}*x'.format(c)
                    elif c < 0:
                        line = line + ' - {0:0.4f}*x'.format(abs(c))
                    else:
                        line = line + ' + {0:0.4f}*x'.format(c)
                else:
                    if line == '':
                        line = line + '{0:0.4f}*x<sup>{1}</sup>'.format(c,i)
                    elif c < 0:
                        line = line + ' - {0:0.4f}*x<sup>{1}</sup>'.format(abs(c),i)
                    else:
                        line = line + ' + {0:0.4f}*x<sup>{1}</sup>'.format(c,i)
                i+=1

        output = output + '<tr>\n<td><u>{0:0.4f} < x <= {1:0.4f}:</u></td>\n</tr>\n<tr>\n<td><b>{2}</b></td>\n</tr>\n'.format(func[1][0],func[1][1],line)
    output = output + '</table>\n'
    return output

def poly_eval(c_list,x):

    i = 0
    res=0
    if all(c == 0 for c in c_list):
        pass
    else:
        for c in c_list:
            res = res + c*math.pow(x,i)
            i+=1

    return res

class no_load:
    def __init__(self, L):
        self.p = 0
        self.rl = 0
        self.rr = 0
        self.L = L

        self.kind = 'NL'

        self.x_graph = [0]
        self.y_graph = [0]

    def chart_load(self,x_scale=0, y_scale=0, arrows=0):
        x = [0]
        y = [0]
        return x,y

    def piece_functions(self):
        '''
        Returns the general piecwise function in the form of two lists
        # list1 is the polynomial coeficients of order [c0,c1x,c2x^2,...,cnx^n]
        # where the list values will only by the cn's*
        # list 2 will be the range over which the function piece applies
        # 0 <= a would be [0,a] **note it will be assumed the the eqality is <= not <
        # rerturned lists will be [[[list11],[list21]],....,[[list1n],[list2n]]
        # where n is the total number of functions to capture the range from
        # 0 to the full span, L of the beam
        '''

        v = [[[0],[0,self.L]]]
        m = [[[0],[0,self.L]]]
        eis = [[[0],[0,self.L]]]
        eid = [[[0],[0,self.L]]]

        vs = PieceFunctionString(v)
        ms = PieceFunctionString(m)
        eiss = PieceFunctionString(eis)
        eids = PieceFunctionString(eid)

        return [v,m,eis,eid],[vs,ms,eiss,eids]

    def fef(self):
        # Fixed End Forces
        RL = 0
        RR = 0
        ML = 0
        MR = 0

        return [RL,ML,RR,MR]

    def v(self,x):
        iters = len(x)
        v=zeros(iters)
        return v

    def m(self,x):
        iters = len(x)
        m=zeros(iters)
        return m

    def eis(self,x):
        iters = len(x)
        eis=zeros(iters)
        return eis

    def eid(self,x):
        iters = len(x)
        eid=zeros(iters)
        return eid

    def vx(self,x):
        v = 0
        return v

    def mx(self,x):
        m = 0
        return m

    def eisx(self,x):
        eisx = 0
        return eisx

    def eidx(self,x):
        eid = 0
        return eid

class pl:
    def __init__(self, p, a, L):

        self.p = float(p)
        self.a = float(a)
        self.L = float(L)
        self.b = self.L - self.a

        self.kind = 'Point'

        self.error = ''

        if self.a > self.L:
            self.error = 'Error a > l'
            self.error = 'Error a > l'

        self.rl = (self.p*self.b)/self.L
        self.rr = (self.p*self.a)/self.L
        self.c4 = ((-1*self.rl * self.a ** 3) / 3) - ((self.rr * self.a ** 3) / 3) + ((self.rr * self.L * self.a ** 2) / 2)
        self.c2 = (-1 / self.L) * ((self.c4) + ((self.rr * self.L ** 3) / 3))
        self.c1 = ((-1*self.rr * self.a ** 2) / 2) - ((self.rl * self.a ** 2) / 2) + (self.rr * self.L * self.a) + self.c2

        arrow_height = self.p/6.0
        #30 degree arrow
        arrow_plus= self.a+(arrow_height*math.tan(math.radians(30)))
        arrow_minus= self.a-(arrow_height*math.tan(math.radians(30)))

        self.x_graph=[arrow_minus,self.a,arrow_plus,self.a,self.a]
        self.y_graph=[arrow_height,0,arrow_height,0,self.p]

    def chart_load(self, x_scale=0, y_scale=0, arrows=0):
        if arrows == 1:
            arrow_height = (self.p/6.0)
            #30 degree arrow
            arrow_plus= (self.a+(arrow_height*math.tan(math.radians(30))))
            arrow_minus= (self.a-(arrow_height*math.tan(math.radians(30))))

            x=[arrow_minus,self.a,arrow_plus,self.a,self.a]
            x = [i*x_scale for i in x]
            y=[arrow_height,0,arrow_height,0,self.p]
            y = [j*y_scale for j in y]
        else:
            x = [self.a*x_scale, self.a*x_scale]
            y = [0,self.p*y_scale]

        return x,y

    def piece_functions(self):
        '''
        Returns the general piecwise function in the form of two lists
        # list1 is the polynomial coeficients of order [c0,c1x,c2x^2,...,cnx^n]
        # where the list values will only by the cn's*
        # list 2 will be the range over which the function piece applies
        # 0 <= a would be [0,a] **note it will be assumed the the eqality is <= not <
        # rerturned lists will be [[[list11],[list21]],....,[[list1n],[list2n]]
        # where n is the total number of functions to capture the range from
        # 0 to the full span, L of the beam
        '''

        if self.a == 0 or self.a == self.L:
            v = [[[0],[0,self.L]]]
            m = [[[0],[0,self.L]]]
            eis = [[[0],[0,self.L]]]
            eid = [[[0],[0,self.L]]]
        else:
            v = [[[self.rl],[0,self.a]],[[-1*self.rr],[self.a,self.L]]]

            m = [[[0,self.rl],[0,self.a]],[[(self.rr * self.L),(-1 * self.rr)],[self.a,self.L]]]

            eis = [[[self.c1,0,self.rl/2.0],[0,self.a]],[[self.c2,(self.rr * self.L),-1.0*self.rr/2.0],[self.a,self.L]]]

            eid = [[[0,self.c1,0,self.rl/6.0],[0,self.a]],[[self.c4, self.c2, self.rr*self.L*0.5,-1*self.rr/6.0],[self.a,self.L]]]

        vs = PieceFunctionString(v)
        ms = PieceFunctionString(m)
        eiss = PieceFunctionString(eis)
        eids = PieceFunctionString(eid)

        return [v,m,eis,eid],[vs,ms,eiss,eids]

    def fef(self):
        # Fixed End Forces
        RL = ((self.p*self.b*self.b) / (self.L*self.L*self.L))*((3*self.a)+self.b)
        RR = ((self.p*self.a*self.a) / (self.L*self.L*self.L))*(self.a+(3*self.b))
        ML = -1*(self.p*self.a*self.b*self.b) / (self.L*self.L)
        MR = (self.p*self.a*self.a*self.b) / (self.L*self.L)

        return [RL,ML,RR,MR]

    def v(self,x):
        iters = len(x)
        v=zeros(iters)

        for i in range(0,iters):
            if x[i] <= self.a:
                if x[i] == 0 and self.a == 0:
                    v[i] = 0
                else:
                    v[i] = self.rl
            else:
                v[i] = -1 * self.rr
        return v

    def m(self,x):
        iters = len(x)
        m=zeros(iters)

        for i in range(0,iters):
            if x[i] <= self.a:
                m[i] = self.rl * x[i]
            else:
                m[i] = (-1 * self.rr * x[i]) + (self.rr * self.L)
        return m

    def eis(self,x):
        iters = len(x)
        eis=zeros(iters)

        for i in range(0,iters):
            if x[i] <= self.a:
                eis[i] = ((self.rl * x[i] ** 2)  / 2) + self.c1
            else:
                eis[i] = ((-1.0 * self.rr * x[i] ** 2)/2.0) + (self.rr * self.L * x[i]) + self.c2
        return eis

    def eid(self,x):
        iters = len(x)
        eid=zeros(iters)

        for i in range(0,iters):
            if x[i] <= self.a:
                eid[i] = ((self.rl * x[i] ** 3) / 6) + (self.c1 * x[i])
            else:
                eid[i] = ((-1*self.rr * x[i] ** 3) / 6) + ((self.rr * self.L * x[i] ** 2) / 2) + (self.c2 * x[i]) + self.c4
        return eid

    def vx(self,x):
        x = float(x)
        if x <= self.a:
            if x==0 and self.a==0:
                v = 0
            else:
                v = self.rl
        else:
            v = -1 * self.rr
        return v

    def mx(self,x):
        x = float(x)
        if x <= self.a:
            m = self.rl * x
        else:
            m = (-1 * self.rr * x) + (self.rr * self.L)
        return m

    def eisx(self,x):
        x = float(x)
        if x <= self.a:
            eisx = ((self.rl * x ** 2)  / 2) + self.c1
        else:
            eisx = ((-1.0 * self.rr * x ** 2)/2.0) + (self.rr * self.L * x) + self.c2
        return eisx

    def eidx(self,x):
        x = float(x)
        if x <= self.a:
            eid = ((self.rl * x ** 3) / 6) + (self.c1 * x)
        else:
            eid = ((-1*self.rr * x ** 3) / 6) + ((self.rr * self.L * x ** 2) / 2) + (self.c2 * x) + self.c4
        return eid

class point_moment:
    def __init__(self, ma, a, L):
        self.ma = float(ma)
        self.a = float(a)
        self.L = float(L)

        self.kind = 'Moment'

        self.error = ''

        if a > self.L:
            self.error = 'Error a > L'

        self.rr = self.ma/self.L
        self.rl = -1.0*self.rr

        self.c2 = (-1.0/self.L) * ((self.ma*self.a**2) - (0.5*self.ma*self.a**2) + (self.rl * (self.L**3/6.0)) + (0.5*self.ma*self.L**2))
        self.c1 = ma*a + self.c2
        self.c3 = 0
        self.c4 = ((-1.0*self.rl*self.L**3)/6.0) - (0.5*self.ma*self.L**2) - (self.c2*self.L)

        r = (self.ma/2.0)
        arrow_height = r/6.0
        #30 degree arrow
        arrow_minus= (arrow_height*math.tan(math.radians(30)))

        if self.ma <0:
            self.x_graph = [self.a,self.a,self.a]
            self.y_graph = [r,0,-r]
            x=0
            y=0
            for a in range(-90, 181):
                x = self.a+(r*math.cos(math.radians(a)))
                y = 0+(r*math.sin(math.radians(a)))
                self.x_graph.append(x)
                self.y_graph.append(y)

            self.x_graph.append(x-arrow_minus)
            self.y_graph.append(y+arrow_height)
            self.x_graph.append(x)
            self.y_graph.append(y)
            self.x_graph.append(x+arrow_minus)
            self.y_graph.append(y+arrow_height)
        else:
            self.x_graph = [self.a-r,self.a,self.a+r, self.a+r-arrow_minus,self.a+r,self.a+r+arrow_minus,self.a+r]
            self.y_graph = [0,0,0,arrow_height,0,arrow_height,0]
            x=0
            y=0
            for a in range(0,271):
                x = self.a+(r*math.cos(math.radians(a)))
                y = 0+(r*math.sin(math.radians(a)))
                self.x_graph.append(x)
                self.y_graph.append(y)

    def chart_load(self, x_scale=0, y_scale=0, arrows=0):
        x=[]
        y=[]
        r = (self.ma/2.0)

        if arrows == 1:
            arrow_height = r/6.0
            #30 degree arrow
            arrow_minus= (arrow_height*math.tan(math.radians(30)))

            if self.ma <0:
                x = [self.a,self.a,self.a]
                y = [r,0,-r]
                xi=0
                yi=0
                for a in range(-90, 181):
                    xi = (self.a)+((r*math.cos(math.radians(a))))
                    yi = 0+((r*math.sin(math.radians(a))))
                    x.append(xi)
                    y.append(yi)

                x.append(xi-arrow_minus)
                y.append(yi+arrow_height)
                x.append(xi)
                y.append(yi)
                x.append(xi+arrow_minus)
                y.append(yi+arrow_height)
            else:
                x = [self.a-r,self.a,self.a+r, self.a+r-arrow_minus,self.a+r,self.a+r+arrow_minus,self.a+r]
                y = [0,0,0,arrow_height,0,arrow_height,0]

                xi=0
                yi=0
                for a in range(0,271):
                    xi = self.a+(r*math.cos(math.radians(a)))
                    yi = 0+(r*math.sin(math.radians(a)))
                    x.append(xi)
                    y.append(yi)
        else:
            if self.ma <0:
                x = [self.a,self.a,self.a]
                y = [r,0,-r]
                xi=0
                yi=0
                for a in range(-90, 181):
                    xi = self.a+(r*math.cos(math.radians(a)))
                    yi = 0+(r*math.sin(math.radians(a)))
                    x.append(xi)
                    y.append(yi)
            else:
                x = [self.a-r,self.a,self.a+r]
                y = [0,r,0]
                xi=0
                yi=0
                for a in range(0,271):
                    xi = self.a+(r*math.cos(math.radians(a)))
                    yi = 0+(r*math.sin(math.radians(a)))
                    x.append(xi)
                    y.append(yi)

        x = [i*x_scale for i in x]
        y = [j*y_scale for j in y]
        return x,y

    def piece_functions(self):
        '''
        Returns the general piecwise function in the form of two lists
        # list1 is the polynomial coeficients of order [c0,c1x,c2x^2,...,cnx^n]
        # where the list values will only by the cn's*
        # list 2 will be the range over which the function piece applies
        # 0 <= a would be [0,a] **note it will be assumed the the eqality is <= not <
        # rerturned lists will be [[[list11],[list21]],....,[[list1n],[list2n]]
        # where n is the total number of functions to capture the range from
        # 0 to the full span, L of the beam
        '''

        v = [[[self.rl],[0,self.L]]]

        if self.a == 0:
            m = [[[self.ma,self.rl],[0,self.L]]]
        elif self.a == self.L:
            m = [[[0,self.rl],[0,self.L]]]
        else:
            m = [[[0,self.rl],[0,self.a]],[[self.ma,self.rl],[self.a,self.L]]]

        eis = [[[self.c1,0,0.5*self.rl],[0,self.a]],[[self.c2,self.ma,0.5*self.rl],[self.a,self.L]]]

        eid = [[[self.c3, self.c1,0,((1/6.0)*self.rl)],[0,self.a]],[[self.c4,self.c2,0.5*self.ma,(1/6.0)*self.rl],[self.a,self.L]]]

        vs = PieceFunctionString(v)
        ms = PieceFunctionString(m)
        eiss = PieceFunctionString(eis)
        eids = PieceFunctionString(eid)

        return [v,m,eis,eid],[vs,ms,eiss,eids]

    def fef(self):
        # Fixed End Forces
        RL = ((-6.0*self.ma*self.a) / (self.L*self.L*self.L)) * (self.L-self.a)
        RR = -1.0*RL
        ML = ((-1.0*self.ma) / (self.L*self.L))*((self.L*self.L)-(4*self.L*self.a)+(3*self.a*self.a))
        MR = -1.0*(self.ma / (self.L*self.L))*((3*self.a*self.a)-(2*self.a*self.L))

        return [RL,ML,RR,MR]

    def v(self,x):
        iters = len(x)
        v=zeros(iters)

        for i in range(0,iters):
            v[i] = self.rl

        return v

    def m(self,x):
        iters = len(x)
        m=zeros(iters)

        for i in range(0,iters):
            if x[i] <= self.a:
                if x[i] == 0 and self.a == 0:
                    m[i] = self.ma
                elif x[i] == self.L and self.a == self.L:
                    m[i] = -1.0*self.ma
                else:
                    m[i] = self.rl * x[i]
            else:
                m[i] = (self.rl * x[i]) + self.ma
        return m

    def eis(self,x):
        iters = len(x)
        eis=zeros(iters)

        for i in range(0,iters):
            if x[i] <= self.a:
                eis[i] = (0.5*self.rl*x[i]**2) + self.c1
            else:
                eis[i] = (0.5*self.rl*x[i]**2) + (self.ma*x[i]) + self.c2
        return eis

    def eid(self,x):
        iters = len(x)
        eid=zeros(iters)

        for i in range(0,iters):
            if x[i] <= self.a:
                eid[i] = ((1/6.0)*self.rl*x[i]**3) + (self.c1*x[i]) + self.c3
            else:
                eid[i] = (1/6.0)*self.rl*x[i]**3 + (0.5*self.ma*x[i]**2) + (self.c2*x[i]) + self.c4
        return eid

    def vx(self,x):
        x = float(x)
        v = self.rl

        return v

    def mx(self,x):
        x = float(x)
        if x <= self.a:
            if x == 0 and self.a == 0:
                m = self.ma
            elif x == self.L and self.a == self.L:
                m = -1.0*self.ma
            else:
                m = self.rl * x
        else:
            m = (self.rl * x) + self.ma
        return m

    def eisx(self,x):
        x = float(x)
        if x <= self.a:
            eis = (0.5*self.rl*x**2) + self.c1
        else:
            eis = (0.5*self.rl*x**2) + (self.ma*x) + self.c2
        return eis

    def eidx(self,x):
        x = float(x)
        if x <= self.a:
            eid = ((1/6.0)*self.rl*x**3) + (self.c1*x) + self.c3
        else:
            eid = (1/6.0)*self.rl*x**3 + (0.5*self.ma*x**2) + (self.c2*x) + self.c4
        return eid

class udl:
    def __init__(self, w1, a, b, L):

        self.w1 = float(w1)
        self.a = float(a)
        self.L = float(L)
        self.b = float(b)
        self.c = b-a

        self.kind = 'UDL'

        self.error = ''

        if self.a > self.b:
            self.error = 'Error a > b'
            self.error = 'Error a > b'
        elif self.a > self.L:
            self.error = 'Error a > l'
            self.error = 'Error a > l'
        elif self.b > self.L:
            self.error = 'Error b > l'
            self.error = 'Error b > l'
        else:
            pass

        self.rl = (self.w1 * self.c) - (((self.w1 * self.c) * (self.a + (self.c / 2))) / self.L)
        self.rr = (((self.w1 * self.c) * (self.a + (self.c / 2))) / self.L)
        self.c1 = 0
        self.c2 = ((-1 * self.w1 * self.a ** 2) / 2)
        self.c3 = self.rr * self.L
        self.c7 = 0
        self.c8 = ((-1 * self.c1 * self.a ** 2) / 2) + ((self.c2 * self.a ** 2) / 2) + ((5 * self.w1 * self.a ** 4) / 24) + self.c7
        self.c9 = ((-1 * self.rl * self.b ** 3) / 3) - ((self.rr * self.b ** 3) / 3) + ((self.w1 * self.b ** 4) / 8) - ((self.w1 * self.a * self.b ** 3) / 3) - ((self.c2 * self.b ** 2) / 2) + ((self.c3 * self.b ** 2) / 2) + self.c8
        self.c6 = ((self.rr * self.L ** 2) / 6) - ((self.c3 * self.L) / 2) - (self.c9 / self.L)
        self.c5 = ((-1 * self.rl * self.b ** 2) / 2) + ((self.w1 * self.b ** 3) / 6) - ((self.w1 * self.a * self.b ** 2) / 2) - ((self.rr * self.b ** 2) / 2) + (self.c3 * self.b) - (self.c2 * self.b) + self.c6
        self.c4 = ((self.w1 * self.a ** 3) / 3) + (self.c2 * self.a) + self.c5 - (self.c1 * self.a)

        arrow_height = self.w1/12.0
        #30 degree arrow
        arrow_plus_start= self.a+(arrow_height*math.tan(math.radians(30)))
        arrow_minus_start= self.a-(arrow_height*math.tan(math.radians(30)))
        arrow_plus_end= self.b+(arrow_height*math.tan(math.radians(30)))
        arrow_minus_end= self.b-(arrow_height*math.tan(math.radians(30)))

        self.x_graph=[arrow_minus_start,self.a,arrow_plus_start,self.a,self.a,self.b,self.b,arrow_minus_end,self.b,arrow_plus_end]
        self.y_graph=[arrow_height,0,arrow_height,0,self.w1,self.w1,0,arrow_height,0,arrow_height]

    def chart_load(self, x_scale=0, y_scale=0, arrows=0):
        x=[]
        y=[]
        if arrows == 1:
            arrow_height = self.w1/6.0
            #30 degree arrow
            arrow_plus_start= self.a+(arrow_height*math.tan(math.radians(30)))
            arrow_minus_start= self.a-(arrow_height*math.tan(math.radians(30)))
            arrow_plus_end= self.b+(arrow_height*math.tan(math.radians(30)))
            arrow_minus_end= self.b-(arrow_height*math.tan(math.radians(30)))

            x=[arrow_minus_start,self.a,arrow_plus_start,self.a,self.a,self.b,self.b,arrow_minus_end,self.b,arrow_plus_end]
            x = [i*x_scale for i in x]
            y=[arrow_height,0,arrow_height,0,self.w1,self.w1,0,arrow_height,0,arrow_height]
            y = [j*y_scale for j in y]
        else:
            x=[self.a,self.a,self.b,self.b]
            x = [i*x_scale for i in x]
            y=[0,self.w1,self.w1,0]
            y = [j*y_scale for j in y]

        return x,y
    def piece_functions(self):
        '''
        Returns the general piecwise function in the form of two lists
        # list1 is the polynomial coeficients of order [c0,c1x,c2x^2,...,cnx^n]
        # where the list values will only by the cn's*
        # list 2 will be the range over which the function piece applies
        # 0 <= a would be [0,a] **note it will be assumed the the eqality is <= not <
        # rerturned lists will be [[[list11],[list21]],....,[[list1n],[list2n]]
        # where n is the total number of functions to capture the range from
        # 0 to the full span, L of the beam
        '''

        v = [[[self.rl],[0,self.a]],[[(self.rl+self.w1*self.a),-1.0*self.w1],[self.a,self.b]],[[-1*self.rr],[self.b,self.L]]]

        m = [[[self.c1,self.rl],[0,self.a]],[[self.c2,self.rl+(self.w1*self.a),-0.5*self.w1],[self.a,self.b]],[[self.c3,-1.0*self.rr],[self.b,self.L]]]

        eis = [[[self.c4,self.c1,0.5*self.rl],[0,self.a]],[[self.c5,self.c2,0.5*(self.rl+(self.w1*self.a)),(-1/6.0)*self.w1],[self.a,self.b]],[[self.c6,self.c3,-0.5*self.rr],[self.b,self.L]]]

        eid = [[[self.c7,self.c4,0.5*self.c1,1/6.0*self.rl],[0,self.a]],[[self.c8, self.c5, 0.5*self.c2,(1/6.0)*(self.rl+(self.w1*self.a)),-1.0*(self.w1 / 24.0)],[self.a,self.b]],[[self.c9,self.c6,0.5*self.c3,((-1.0 * self.rr) / 6.0)],[self.b,self.L]]]

        vs = PieceFunctionString(v)
        ms = PieceFunctionString(m)
        eiss = PieceFunctionString(eis)
        eids = PieceFunctionString(eid)

        return [v,m,eis,eid],[vs,ms,eiss,eids]

    def v(self,x):
        iters = len(x)
        v=zeros(iters)

        for i in range(0,iters):
            if x[i] <= self.a:
                v[i] = self.rl
            elif x[i]<=self.b:
                v[i] = self.rl - (self.w1 * (x[i] - self.a))
            else:
                v[i] = -1 * self.rr
        return v

    def m(self,x):
        iters = len(x)
        m=zeros(iters)

        for i in range(0,iters):
            if x[i] <= self.a:
                m[i] = (self.rl * x[i]) + self.c1
            elif x[i] <= self.b:
                m[i] = (self.rl * x[i]) - ((self.w1 * x[i] ** 2) / 2) + (self.w1 * self.a * x[i]) + self.c2
            else:
                m[i] = (-1 * self.rr * x[i]) + self.c3
        return m

    def eis(self,x):
        iters = len(x)
        eis=zeros(iters)

        for i in range(0,iters):
            if x[i] <= self.a:
                eis[i] = ((self.rl * x[i] ** 2) / 2.0) + (self.c1 * x[i]) + self.c4
            elif x[i] <= self.b:
                eis[i] = ((self.rl * x[i] **2) / 2.0) - ((self.w1 * x[i] ** 3) / 6.0) + ((self.w1 * self.a * x[i] **2) / 2.0) + (self.c2 * x[i]) + self.c5
            else:
                eis[i] = ((-1.0 * self.rr * x[i] ** 2) / 2.0) + (self.c3 * x[i]) + self.c6
        return eis

    def eid(self,x):
        iters = len(x)
        eid=zeros(iters)

        for i in range(0,iters):
            if x[i] <= self.a:
                eid[i] = ((self.rl * x[i] ** 3) / 6) + ((self.c1 * x[i] ** 2) / 2) + (self.c4 * x[i]) + self.c7
            elif x[i]<=self.b:
                eid[i] = ((self.rl * x[i] ** 3) / 6) - ((self.w1 * x[i] ** 4) / 24) + ((self.w1 * self.a * x[i] ** 3) / 6) + ((self.c2 * x[i] ** 2) / 2) + (self.c5 * x[i]) + self.c8
            else:
                eid[i] = ((-1 * self.rr * x[i] ** 3) / 6) + ((self.c3 * x[i] ** 2) / 2) + (self.c6 * x[i]) + self.c9
        return eid

    def vx(self,x):
        x = float(x)
        if x <= self.a:
            v = self.rl
        elif x<=self.b:
            v = self.rl - (self.w1 * (x - self.a))
        else:
            v = -1 * self.rr
        return v

    def mx(self,x):
        x = float(x)
        if x <= self.a:
            m = (self.rl * x) + self.c1
        elif x <= self.b:
            m = (self.rl * x) - ((self.w1 * x ** 2) / 2) + (self.w1 * self.a * x) + self.c2
        else:
            m = (-1 * self.rr * x) + self.c3
        return m

    def eisx(self,x):
        x = float(x)
        if x <= self.a:
            eis = ((self.rl * x ** 2) / 2.0) + (self.c1 * x) + self.c4
        elif x <= self.b:
            eis = ((self.rl * x **2) / 2.0) - ((self.w1 * x ** 3) / 6.0) + ((self.w1 * self.a * x **2) / 2.0) + (self.c2 * x) + self.c5
        else:
            eis = ((-1.0 * self.rr * x ** 2) / 2.0) + (self.c3 * x) + self.c6
        return eis

    def eidx(self,x):
        x = float(x)
        if x <= self.a:
            eid = ((self.rl * x ** 3) / 6) + ((self.c1 * x ** 2) / 2) + (self.c4 * x) + self.c7
        elif x<=self.b:
            eid = ((self.rl * x ** 3) / 6) - ((self.w1 * x ** 4) / 24) + ((self.w1 * self.a * x ** 3) / 6) + ((self.c2 * x ** 2) / 2) + (self.c5 * x) + self.c8
        else:
            eid = ((-1 * self.rr * x ** 3) / 6) + ((self.c3 * x ** 2) / 2) + (self.c6 * x) + self.c9
        return eid

    def fef(self):
        eis0 = self.eisx(0)
        eisL = self.eisx(self.L)

        s = np.array([[-1.0*eis0],[-1.0*eisL]])

        ems = np.array([[-1.0*self.L/3.0 , self.L/6.0],[self.L/6.0 , -1.0*self.L/3.0]])

        fem = np.linalg.solve(ems,s)

        mo = point_moment(fem[0][0],0,self.L)
        ml = point_moment(fem[1][0],self.L,self.L)

        RL = self.rl+mo.rl+ml.rl
        RR = self.rr+mo.rr+ml.rr
        ML = fem[0][0]
        MR = fem[1][0]

        return [RL,ML,RR,MR]

class trap:
    def __init__(self, w1, w2, a, b, L):

        self.w1 = float(w1)
        self.w2 = float(w2)
        self.a = float(a)
        self.L = float(L)
        self.b = float(b)
        self.c = self.b-self.a

        self.kind = 'TRAP'


        self.error = ''

        if self.a > self.b:
            self.error = 'Error a > b'
            self.error = 'Error a > b'
        elif self.a > self.L:
            self.error = 'Error a > l'
            self.error = 'Error a > l'
        elif self.b > self.L:
            self.error = 'Error b > l'
            self.error = 'Error b > l'
        elif sign(self.w1) != sign(self.w2) and self.w1 !=0 and self.w2 !=0:
            self.error = 'Error w1 and w2 change direction'
            self.error = 'Error w1 and w2 change direction'
        else:
            pass

        self.s = (self.w2 -self.w1)/self.c
        self.xbar = (self.c * ((2 * self.w2) + self.w1)) / (3 * (self.w2 + self.w1))
        self.W = self.c * ((self.w1 + self.w2) / 2)
        self.rr = (self.W * (self.a + self.xbar)) / self.L
        self.rl = self.W - self.rr
        self.c1 = 0
        self.c2 = self.c1 + ((self.a ** 3 * self.s) / 6) + ((self.a ** 2 * (self.w1 - (self.s * self.a))) / 2) + ((((self.s * self.a) - (2 * self.w1)) * self.a ** 2) / 2)
        self.c3 = self.rr * self.L
        self.c7 = 0
        self.c8 = ((-1 * self.c1 * self.a ** 2) / 2) - ((self.a ** 5 * self.s) / 30) - ((self.a ** 4 * (self.w1 - (self.s * self.a))) / 8) - ((((self.s * self.a) - (2 * self.w1)) * self.a ** 4) / 6) + ((self.c2 * self.a ** 2) / 2) + self.c7
        self.c9 = ((-1 * self.rl * self.b ** 3) / 3) + ((self.b ** 5 * self.s) / 30) + ((self.b ** 4 * (self.w1 - (self.s * self.a))) / 8) + ((((self.s * self.a) - (2 * self.w1)) * self.a * self.b ** 3) / 6) - ((self.c2 * self.b ** 2) / 2) + self.c8 - ((self.rr * self.b ** 3) / 3) + ((self.c3 * self.b ** 2) / 2)
        self.c6 = (((self.rr * self.L ** 3) / 6) - ((self.c3 * self.L ** 2) / 2) - self.c9) / self.L
        self.c5 = ((-1 * self.rr * self.b ** 2) / 2) + (self.c3 * self.b) + self.c6 - ((self.rl * self.b ** 2) / 2) + ((self.b ** 4 * self.s) / 24) + ((self.b ** 3 * (self.w1 - (self.s * self.a))) / 6) + ((((self.s * self.a) - (2 * self.w1)) * self.a * self.b ** 2) / 4) - (self.c2 * self.b)
        self.c4 = ((-1 * self.a ** 4 * self.s) / 24) - ((self.a ** 3 * (self.w1 - (self.s * self.a))) / 6) - ((((self.s * self.a) - (2 * self.w1)) * self.a ** 3) / 4) + (self.c2 * self.a) + self.c5 - (self.c1 * self.a)

        arrow_height = self.w1/6.0
        arrow_height2 = self.w2/6.0
        #30 degree arrow
        arrow_plus_start= self.a+(arrow_height*math.tan(math.radians(30)))
        arrow_minus_start= self.a-(arrow_height*math.tan(math.radians(30)))
        arrow_plus_end= self.b+(arrow_height2*math.tan(math.radians(30)))
        arrow_minus_end= self.b-(arrow_height2*math.tan(math.radians(30)))

        self.x_graph=[arrow_minus_start,self.a,arrow_plus_start,self.a,self.a,self.b,self.b,arrow_minus_end,self.b,arrow_plus_end]
        self.y_graph=[arrow_height,0,arrow_height,0,self.w1,self.w2,0,arrow_height2,0,arrow_height2]

    def chart_load(self, x_scale=0, y_scale=0, arrows=0):
        x=[]
        y=[]
        if arrows == 1:
            arrow_height = self.w1/6.0
            arrow_height2 = self.w2/6.0
            #30 degree arrow
            arrow_plus_start= self.a+(arrow_height*math.tan(math.radians(30)))
            arrow_minus_start= self.a-(arrow_height*math.tan(math.radians(30)))
            arrow_plus_end= self.b+(arrow_height2*math.tan(math.radians(30)))
            arrow_minus_end= self.b-(arrow_height2*math.tan(math.radians(30)))

            x=[arrow_minus_start,self.a,arrow_plus_start,self.a,self.a,self.b,self.b,arrow_minus_end,self.b,arrow_plus_end]
            x = [i*x_scale for i in x]
            y=[arrow_height,0,arrow_height,0,self.w1,self.w2,0,arrow_height2,0,arrow_height2]
            y = [j*y_scale for j in y]

        else:
            x=[self.a,self.a,self.b,self.b]
            x = [i*x_scale for i in x]
            y=[0,self.w1,self.w2,0]
            y = [j*y_scale for j in y]

        return x,y

    def piece_functions(self):
        '''
        Returns the general piecwise function in the form of two lists
        # list1 is the polynomial coeficients of order [c0,c1x,c2x^2,...,cnx^n]
        # where the list values will only by the cn's*
        # list 2 will be the range over which the function piece applies
        # 0 <= a would be [0,a] **note it will be assumed the the eqality is <= not <
        # rerturned lists will be [[[list11],[list21]],....,[[list1n],[list2n]]
        # where n is the total number of functions to capture the range from
        # 0 to the full span, L of the beam
        '''

        v = [[[self.rl],[0,self.a]],[[self.rl- ((((self.s * self.a) - (2 * self.w1)) * self.a) / 2),-1.0*((self.w1 - (self.s * self.a))),-1.0*(self.s/ 2)],[self.a,self.b]],[[-1.0*self.rr],[self.b,self.L]]]

        m = [[[self.c1,self.rl],[0,self.a]],[[self.c2,self.rl - ((((self.s * self.a) - (2 * self.w1)) * self.a) / 2.0),-1.0*((self.w1 - (self.s * self.a)) / 2.0),-1.0*((self.s) / 6.0)],[self.a,self.b]],[[self.c3,-1.0*self.rr],[self.b,self.L]]]

        eis = [[[self.c4,self.c1,(self.rl / 2.0)],[0,self.a]],[[self.c5,self.c2,(self.rl/ 2.0) - ((((self.s * self.a) - (2 * self.w1)) * self.a) / 4.0), -1.0*((self.w1 - (self.s * self.a)) / 6.0),-1.0*(self.s / 24.0)],[self.a,self.b]],[[self.c6,self.c3,((-1.0* self.rr) / 2)],[self.b,self.L]]]

        eid =  [[[self.c7,self.c4,(self.c1 / 2.0),(self.rl/ 6.0)],[0,self.a]],[[self.c8,self.c5,self.c2 / 2.0,(self.rl / 6.0) - ((((self.s * self.a) - (2 * self.w1)) * self.a) / 12.0), -1.0*((self.w1 - (self.s * self.a)) / 24),-1.0*(self.s / 120.0)],[self.a,self.b]],[[self.c9,self.c6,(self.c3 / 2.0),((-1.0 * self.rr) / 6.0)],[self.b,self.L]]]

        vs = PieceFunctionString(v)
        ms = PieceFunctionString(m)
        eiss = PieceFunctionString(eis)
        eids = PieceFunctionString(eid)

        return [v,m,eis,eid],[vs,ms,eiss,eids]

    def v(self,x):
        iters = len(x)
        v=zeros(iters)

        for i in range(0,iters):
            if x[i] <= self.a:
                v[i] = self.rl
            elif x[i]<=self.b:
                v[i] = self.rl - ((x[i] ** 2 * self.s) / 2) - (x[i] * (self.w1 - (self.s * self.a))) - ((((self.s * self.a) - (2 * self.w1)) * self.a) / 2)
            else:
                v[i] = -1 * self.rr
        return v

    def m(self,x):
        iters = len(x)
        m=zeros(iters)

        for i in range(0,iters):
            if x[i] <= self.a:
                m[i] = (self.rl * x[i]) + self.c1
            elif x[i] <= self.b:
                m[i] = (self.rl * x[i]) - ((x[i] ** 3 * self.s) / 6) - ((x[i] ** 2 * (self.w1 - (self.s * self.a))) / 2) - ((((self.s * self.a) - (2 * self.w1)) * self.a * x[i]) / 2) + self.c2
            else:
                m[i] = (-1 * self.rr * x[i]) + self.c3
        return m

    def eis(self,x):
        iters = len(x)
        eis=zeros(iters)

        for i in range(0,iters):
            if x[i] <= self.a:
                eis[i] = ((self.rl * x[i] ** 2) / 2) + (self.c1 * x[i]) + self.c4
            elif x[i] <= self.b:
                eis[i] = ((self.rl * x[i] ** 2) / 2) - ((x[i] ** 4 * self.s) / 24) - ((x[i] ** 3 * (self.w1 - (self.s * self.a))) / 6) - ((((self.s * self.a) - (2 * self.w1)) * self.a * x[i] ** 2) / 4) + (self.c2 * x[i]) + self.c5
            else:
                eis[i] = ((-1 * self.rr * x[i] ** 2) / 2) + (self.c3 * x[i]) + self.c6
        return eis

    def eid(self,x):
        iters = len(x)
        eid=zeros(iters)

        for i in range(0,iters):
            if x[i] <= self.a:
                eid[i] = ((self.rl * x[i] ** 3) / 6) + ((self.c1 * x[i] ** 2) / 2) + (self.c4 * x[i]) + self.c7
            elif x[i]<=self.b:
                eid[i] = ((self.rl * x[i] ** 3) / 6) - ((x[i] ** 5 * self.s) / 120) - ((x[i] ** 4 * (self.w1 - (self.s * self.a))) / 24) - ((((self.s * self.a) - (2 * self.w1)) * self.a * x[i] ** 3) / 12) + ((self.c2 * x[i] ** 2) / 2) + (self.c5 * x[i]) + self.c8
            else:
                eid[i] = ((-1 * self.rr * x[i] ** 3) / 6) + ((self.c3 * x[i] ** 2) / 2) + (self.c6 * x[i]) + self.c9
        return eid

    def vx(self,x):
        x = float(x)
        if x <= self.a:
            v = self.rl
        elif x<=self.b:
            v = self.rl - ((x ** 2 * self.s) / 2) - (x * (self.w1 - (self.s * self.a))) - ((((self.s * self.a) - (2 * self.w1)) * self.a) / 2)
        else:
            v = -1 * self.rr
        return v

    def mx(self,x):
        x = float(x)
        if x <= self.a:
            m = (self.rl * x) + self.c1
        elif x <= self.b:
            m = (self.rl * x) - ((x ** 3 * self.s) / 6) - ((x ** 2 * (self.w1 - (self.s * self.a))) / 2) - ((((self.s * self.a) - (2 * self.w1)) * self.a * x) / 2) + self.c2
        else:
            m = (-1 * self.rr * x) + self.c3
        return m

    def eisx(self,x):
        x = float(x)
        if x <= self.a:
            eis = ((self.rl * x ** 2) / 2) + (self.c1 * x) + self.c4
        elif x <= self.b:
            eis = ((self.rl * x ** 2) / 2) - ((x ** 4 * self.s) / 24) - ((x ** 3 * (self.w1 - (self.s * self.a))) / 6) - ((((self.s * self.a) - (2 * self.w1)) * self.a * x ** 2) / 4) + (self.c2 * x) + self.c5
        else:
            eis = ((-1 * self.rr * x ** 2) / 2) + (self.c3 * x) + self.c6
        return eis

    def eidx(self,x):
        x = float(x)
        if x <= self.a:
            eid = ((self.rl * x ** 3) / 6) + ((self.c1 * x ** 2) / 2) + (self.c4 * x) + self.c7
        elif x<=self.b:
            eid = ((self.rl * x ** 3) / 6) - ((x ** 5 * self.s) / 120) - ((x ** 4 * (self.w1 - (self.s * self.a))) / 24) - ((((self.s * self.a) - (2 * self.w1)) * self.a * x ** 3) / 12) + ((self.c2 * x ** 2) / 2) + (self.c5 * x) + self.c8
        else:
            eid = ((-1 * self.rr * x ** 3) / 6) + ((self.c3 * x ** 2) / 2) + (self.c6 * x) + self.c9
        return eid

    def fef(self):
        eis0 = self.eisx(0)
        eisL = self.eisx(self.L)

        s = np.array([[-1.0*eis0],[-1.0*eisL]])

        ems = np.array([[-1.0*self.L/3.0 , self.L/6.0],[self.L/6.0 , -1.0*self.L/3.0]])

        fem = np.linalg.solve(ems,s)

        mo = point_moment(fem[0][0],0,self.L)
        ml = point_moment(fem[1][0],self.L,self.L)

        RL = self.rl+mo.rl+ml.rl
        RR = self.rr+mo.rr+ml.rr
        ML = fem[0][0]
        MR = fem[1][0]

        return [RL,ML,RR,MR]

class end_delta:
    def __init__(self, delta_i, delta_j, L):
        '''
        Important note it is assumed that delta_i and delta_j
        have been divided by E and I. If this is being used
        in combination with other loads make sure consistent
        units are being used
        '''

        self.rl = 0
        self.rr = 0
        self.deltai = delta_i
        self.deltaj = delta_j
        self.L = L

        self.slope = (delta_j - delta_i)/self.L

        self.kind = 'END_DELTA'

        self.x_graph = [0]
        self.y_graph = [0]

    def chart_load(self, x_scale=0, y_scale=0, arrows=0):
        x=[0]
        y=[0]

        return x,y
    def piece_functions(self):
        '''
        Returns the general piecwise function in the form of two lists
        # list1 is the polynomial coeficients of order [c0,c1x,c2x^2,...,cnx^n]
        # where the list values will only by the cn's*
        # list 2 will be the range over which the function piece applies
        # 0 <= a would be [0,a] **note it will be assumed the the eqality is <= not <
        # rerturned lists will be [[[list11],[list21]],....,[[list1n],[list2n]]
        # where n is the total number of functions to capture the range from
        # 0 to the full span, L of the beam
        '''

        v = [[[0],[0,self.L]]]
        m = [[[0],[0,self.L]]]
        eis = [[[self.slope],[0,self.L]]]
        eid = [[[self.deltai,self.slope],[0,self.L]]]

        vs = PieceFunctionString(v)
        ms = PieceFunctionString(m)
        eiss = PieceFunctionString(eis)
        eids = PieceFunctionString(eid)

        return [v,m,eis,eid],[vs,ms,eiss,eids]

    def v(self,x):
        iters = len(x)
        v=zeros(iters)
        return v

    def m(self,x):
        iters = len(x)
        m=zeros(iters)
        return m

    def eis(self,x):
        iters = len(x)
        eis=zeros(iters)
        for i in range(0,iters):
            eis[i] = self.slope
        return eis

    def eid(self,x):
        iters = len(x)
        eid=zeros(iters)
        for i in range(0,iters):
            eid[i] = self.slope*x[i] + self.deltai
        return eid

    def vx(self,x):
        v = 0
        return v

    def mx(self,x):
        m = 0
        return m

    def eisx(self,x):
        eisx = self.slope
        return eisx

    def eidx(self,x):
        eid = self.slope*x + self.deltai
        return eid

    def fef(self):
        eis0 = self.eisx(0)
        eisL = self.eisx(self.L)

        s = np.array([[-1.0*eis0],[-1.0*eisL]])

        ems = np.array([[-1.0*self.L/3.0 , self.L/6.0],[self.L/6.0 , -1.0*self.L/3.0]])

        fem = np.linalg.solve(ems,s)

        mo = point_moment(fem[0][0],0,self.L)
        ml = point_moment(fem[1][0],self.L,self.L)

        RL = self.rl+mo.rl+ml.rl
        RR = self.rr+mo.rr+ml.rr
        ML = fem[0][0]
        MR = fem[1][0]

        return [RL,ML,RR,MR]

class cant_right_nl:
    def __init__(self, slope,L):
        self.slope = slope
        self.L = L
        self.rl = 0
        self.rr = 0
        self.ml = 0

        self.kind = 'NL'

        self.x_graph = [0]
        self.y_graph = [0]

    def chart_load(self, x_scale=0, y_scale=0, arrows=0):
        x=[0]
        y=[0]

        return x,y

    def piece_functions(self):
        '''
        Returns the general piecwise function in the form of two lists
        # list1 is the polynomial coeficients of order [c0,c1x,c2x^2,...,cnx^n]
        # where the list values will only by the cn's*
        # list 2 will be the range over which the function piece applies
        # 0 <= a would be [0,a] **note it will be assumed the the eqality is <= not <
        # rerturned lists will be [[[list11],[list21]],....,[[list1n],[list2n]]
        # where n is the total number of functions to capture the range from
        # 0 to the full span, L of the beam
        '''

        v = [[[0],[0,self.L]]]

        m = [[[0],[0,self.L]]]

        eis = [[[self.slope],[0,self.L]]]

        eid = [[[0, self.slope],[0,self.L]]]

        vs = PieceFunctionString(v)
        ms = PieceFunctionString(m)
        eiss = PieceFunctionString(eis)
        eids = PieceFunctionString(eid)

        return [v,m,eis,eid],[vs,ms,eiss,eids]

    def fef(self):
        # Fixed End Forces
        RL = 0
        RR = 0
        ML = 0
        MR = 0

        return [RL,ML,RR,MR]

    def v(self,x):
        iters = len(x)
        v=zeros(iters)

        return v

    def m(self,x):
        iters = len(x)
        m=zeros(iters)

        return m

    def eis(self,x):
        iters = len(x)
        eis=zeros(iters)
        for i in range(0,iters):
            eis[i] = self.slope

        return eis

    def eid(self,x):
        iters = len(x)
        eid=zeros(iters)
        for i in range(0,iters):
            eid[i] = self.slope * x[i]

        return eid

    def vx(self,x):
        v=0

        return v

    def mx(self,x):
        m=0

        return m

    def eisx(self,x):
        eis = self.slope

        return eis

    def eidx(self,x):
        eid = self.slope * x

        return eid

class cant_right_point:
    def __init__(self, p, a, L, Lb):

        self.p = float(p)
        self.a = float(a)
        self.L = float(L)
        self.Lb = float(Lb)
        self.b = self.L - self.a

        self.kind = 'Point'

        self.error = ''

        if self.a > self.L:
            self.error = 'Error a > l'
            self.error = 'Error a > l'

        self.rl = self.p
        self.rr = 0
        self.ml = -1.0*self.p*self.a

        # 0 length backspan indicates fixed-free beam initialize slope to 0
        if Lb == 0:
            self.backspan = no_load(0)
            self.c1 = 0
        else:
            self.backspan = point_moment(-1.0*self.ml,self.Lb,self.Lb)
            self.c1 = self.backspan.eisx(self.Lb)

        self.c2 = 0
        self.c3 = 0.5*self.rl*self.a**2 + self.ml*self.a + self.c1
        self.c4 = -1.0*self.c3*self.a + (1.0/6.0)*self.rl*self.a**3 + 0.5*self.ml*self.a**2 + self.c1*self.a + self.c2

        arrow_height = self.p/6.0
        #30 degree arrow
        arrow_plus= self.a+(arrow_height*math.tan(math.radians(30)))
        arrow_minus= self.a-(arrow_height*math.tan(math.radians(30)))

        self.x_graph=[arrow_minus,self.a,arrow_plus,self.a,self.a]
        self.y_graph=[arrow_height,0,arrow_height,0,self.p]

    def chart_load(self, x_scale=0, y_scale=0, arrows=0):

        if arrows == 1:
            arrow_height = self.p/6.0
            #30 degree arrow
            arrow_plus= self.a+(arrow_height*math.tan(math.radians(30)))
            arrow_minus= self.a-(arrow_height*math.tan(math.radians(30)))

            x=[arrow_minus,self.a,arrow_plus,self.a,self.a]
            x = [i*x_scale for i in x]
            y=[arrow_height,0,arrow_height,0,self.p]
            y = [j*y_scale for j in y]
        else:
            x=[self.a,self.a]
            x = [i*x_scale for i in x]
            y=[0,self.p]
            y = [j*y_scale for j in y]

        return x,y

    def piece_functions(self):
        '''
        Returns the general piecwise function in the form of two lists
        # list1 is the polynomial coeficients of order [c0,c1x,c2x^2,...,cnx^n]
        # where the list values will only by the cn's*
        # list 2 will be the range over which the function piece applies
        # 0 <= a would be [0,a] **note it will be assumed the the eqality is <= not <
        # rerturned lists will be [[[list11],[list21]],....,[[list1n],[list2n]]
        # where n is the total number of functions to capture the range from
        # 0 to the full span, L of the beam
        '''

        if self.a == 0:
            v = [[[0],[0,self.L]]]
            m = [[[0],[0,self.L]]]
            eis = [[[0],[0,self.L]]]
            eid = [[[0],[0,self.L]]]
        else:
            v = [[[self.p],[0,self.a]],[[0],[self.a,self.L]]]

            m = [[[self.ml,self.rl],[0,self.a]],[[0],[self.a,self.L]]]

            eis = [[[self.c1,self.ml,0.5*self.rl],[0,self.a]],[[self.c3],[self.a,self.L]]]

            eid = [[[self.c2,self.c1,0.5*self.ml,(1.0/6.0)*self.rl],[0,self.a]],[[self.c4, self.c3],[self.a,self.L]]]

        vs = PieceFunctionString(v)
        ms = PieceFunctionString(m)
        eiss = PieceFunctionString(eis)
        eids = PieceFunctionString(eid)

        return [v,m,eis,eid],[vs,ms,eiss,eids]

    def fef(self):
        # Fixed End Forces
        RL = self.rl
        RR = 0
        ML = self.ml
        MR = 0

        return [RL,ML,RR,MR]

    def v(self,x):
        iters = len(x)
        v=zeros(iters)

        for i in range(0,iters):
            if x[i]<=self.a:
                if x[i] == 0 and self.a == 0:
                    v[i] == 0
                else:
                    v[i] = self.p
            else:
                v[i] = 0
        return v

    def m(self,x):
        iters = len(x)
        m=zeros(iters)

        for i in range(0,iters):
            if x[i]<=self.a:
                m[i] = self.rl*x[i] + self.ml
            else:
                m[i] = 0
        return m

    def eis(self,x):
        iters = len(x)
        eis=zeros(iters)

        for i in range(0,iters):
            if x[i]<=self.a:
                eis[i] = 0.5*self.rl*x[i]**2 + self.ml*x[i] + self.c1
            else:
                eis[i] = self.c3
        return eis

    def eid(self,x):
        iters = len(x)
        eid=zeros(iters)

        for i in range(0,iters):
            if x[i]<=self.a:
                eid[i] = (1.0/6.0)*self.rl*x[i]**3 + 0.5*self.ml*x[i]**2 + self.c1*x[i] + self.c2
            else:
                eid[i] = self.c3*x[i] + self.c4
        return eid

    def vx(self,x):
        if x<=self.a:
            if x == 0 and self.a ==0:
                v = 0
            else:
                v = self.p
        else:
            v = 0
        return v

    def mx(self,x):
        if x<=self.a:
            m = self.rl*x + self.ml
        else:
            m = 0
        return m

    def eisx(self,x):
        if x<=self.a:
            eis = 0.5*self.rl*x**2 + self.ml*x + self.c1
        else:
            eis = self.c3
        return eis

    def eidx(self,x):
        if x<=self.a:
            eid = (1.0/6.0)*self.rl*x**3 + 0.5*self.ml*x**2 + self.c1*x + self.c2
        else:
            eid = self.c3*x + self.c4
        return eid

class cant_right_point_moment:
    def __init__(self, ma, a, L, Lb):

        self.ma = float(ma)
        self.a = float(a)
        self.L = float(L)
        self.Lb = float(Lb)
        self.b = self.L - self.a

        self.kind = 'Moment'

        self.error = ''

        if self.a > self.L:
            self.error = 'Error a > l'
            self.error = 'Error a > l'

        self.rl = 0
        self.rr = 0
        self.ml = -1.0*self.ma

        # 0 length backspan indicates fixed-free beam initialize slope to 0
        if Lb == 0:
            self.backspan = no_load(0)
            self.c1 = 0
        else:
            self.backspan = point_moment(-1.0*self.ml,self.Lb,self.Lb)
            self.c1 = self.backspan.eisx(self.Lb)

        self.c2 = 0
        self.c3 = self.ml*self.a + self.c1
        self.c4 = 0.5*self.ml*self.a**2 + self.c1 * self.a + self.c2 - self.c3 * self.a

        r = (self.ma/2.0)
        arrow_height = r/6.0
        #30 degree arrow
        arrow_minus= (arrow_height*math.tan(math.radians(30)))

        if self.ma <0:
            self.x_graph = [self.a,self.a,self.a]
            self.y_graph = [r,0,-r]
            x=0
            y=0
            for a in range(-90, 181):
                x = self.a+(r*math.cos(math.radians(a)))
                y = 0+(r*math.sin(math.radians(a)))
                self.x_graph.append(x)
                self.y_graph.append(y)

            self.x_graph.append(x-arrow_minus)
            self.y_graph.append(y+arrow_height)
            self.x_graph.append(x)
            self.y_graph.append(y)
            self.x_graph.append(x+arrow_minus)
            self.y_graph.append(y+arrow_height)
        else:
            self.x_graph = [self.a-r,self.a,self.a+r, self.a+r-arrow_minus,self.a+r,self.a+r+arrow_minus,self.a+r]
            self.y_graph = [0,0,0,arrow_height,0,arrow_height,0]
            x=0
            y=0
            for a in range(0,271):
                x = self.a+(r*math.cos(math.radians(a)))
                y = 0+(r*math.sin(math.radians(a)))
                self.x_graph.append(x)
                self.y_graph.append(y)

    def chart_load(self, x_scale=0, y_scale=0, arrows=0):
        r = (self.ma/2.0)

        if arrows == 1:
            arrow_height = r/6.0
            #30 degree arrow
            arrow_minus= (arrow_height*math.tan(math.radians(30)))

            if self.ma <0:
                x= [self.a,self.a,self.a]
                y = [r,0,-r]
                xi=0
                yi=0
                for a in range(-90, 181):
                    xi = self.a+(r*math.cos(math.radians(a)))
                    yi = 0+(r*math.sin(math.radians(a)))
                    x.append(xi)
                    y.append(yi)

                x.append(xi-arrow_minus)
                y.append(yi+arrow_height)
                x.append(xi)
                y.append(yi)
                x.append(xi+arrow_minus)
                y.append(yi+arrow_height)
            else:
                x= [self.a-r,self.a,self.a+r, self.a+r-arrow_minus,self.a+r,self.a+r+arrow_minus,self.a+r]
                y = [0,0,0,arrow_height,0,arrow_height,0]
                xi=0
                yi=0
                for a in range(0,271):
                    xi = self.a+(r*math.cos(math.radians(a)))
                    yi = 0+(r*math.sin(math.radians(a)))
                    x.append(xi)
                    y.append(yi)

        x = [i*x_scale for i in x]
        y = [j*y_scale for j in y]

        return x,y

    def piece_functions(self):
        '''
        Returns the general piecwise function in the form of two lists
        # list1 is the polynomial coeficients of order [c0,c1x,c2x^2,...,cnx^n]
        # where the list values will only by the cn's*
        # list 2 will be the range over which the function piece applies
        # 0 <= a would be [0,a] **note it will be assumed the the eqality is <= not <
        # rerturned lists will be [[[list11],[list21]],....,[[list1n],[list2n]]
        # where n is the total number of functions to capture the range from
        # 0 to the full span, L of the beam
        '''

        v = [[[0],[0,self.L]]]

        m = [[[self.ml],[0,self.a]],[[0],[self.a,self.L]]]

        eis = [[[self.c1,self.ml],[0,self.a]],[[self.c3],[self.a,self.L]]]

        eid = [[[self.c2, self.c1,0.5*self.ml],[0,self.a]],[[self.c4,self.c3],[self.a,self.L]]]

        vs = PieceFunctionString(v)
        ms = PieceFunctionString(m)
        eiss = PieceFunctionString(eis)
        eids = PieceFunctionString(eid)

        return [v,m,eis,eid],[vs,ms,eiss,eids]

    def fef(self):
        # Fixed End Forces
        RL = self.rl
        RR = 0
        ML = self.ml
        MR = 0

        return [RL,ML,RR,MR]

    def v(self,x):
        iters = len(x)
        v=zeros(iters)

        for i in range(0,iters):
            if x[i]<=self.a:
                v[i] = 0
            else:
                v[i] = 0
        return v

    def m(self,x):
        iters = len(x)
        m=zeros(iters)

        for i in range(0,iters):
            if x[i]<=self.a:
                m[i] = self.ml
            else:
                m[i] = 0
        return m

    def eis(self,x):
        iters = len(x)
        eis=zeros(iters)

        for i in range(0,iters):
            if x[i]<=self.a:
                eis[i] = self.ml*x[i] + self.c1
            else:
                eis[i] = self.c3
        return eis

    def eid(self,x):
        iters = len(x)
        eid=zeros(iters)

        for i in range(0,iters):
            if x[i]<=self.a:
                eid[i] = 0.5*self.ml*x[i]**2 + self.c1*x[i] + self.c2
            else:
                eid[i] = self.c3*x[i] + self.c4
        return eid

    def vx(self,x):
        if x<=self.a:
            v = 0
        else:
            v = 0
        return v

    def mx(self,x):
        if x<=self.a:
            m = self.ml
        else:
            m = 0
        return m

    def eisx(self,x):
        if x<=self.a:
            eis = self.ml*x + self.c1
        else:
            eis = self.c3
        return eis

    def eidx(self,x):
        if x<=self.a:
            eid = 0.5*self.ml*x**2 + self.c1*x + self.c2
        else:
            eid = self.c3*x + self.c4
        return eid

class cant_right_udl:
    def __init__(self, w1, a, b, L, Lb):

        self.w1 = float(w1)
        self.a = float(a)
        self.L = float(L)
        self.b = float(b)
        self.c = self.b - self.a
        self.w_tot = self.w1*self.c
        self.Lb = float(Lb)

        self.kind = 'UDL'

        self.error = ''

        if self.a > self.b:
            self.error = 'Error a > b'
            self.error = 'Error a > b'
        elif self.a > self.L:
            self.error = 'Error a > l'
            self.error = 'Error a > l'
        elif self.b > self.L:
            self.error = 'Error b > l'
            self.error = 'Error b > l'
        else:
            pass

        self.rl = self.w_tot
        self.rr = 0
        self.ml = -1.0*self.w_tot*(self.b-(self.c/2))

        # 0 length backspan indicates fixed-free beam initialize slope to 0
        if Lb == 0:
            self.backspan = no_load(0)
            self.c1 = 0
        else:
            self.backspan = point_moment(-1.0*self.ml,self.Lb,self.Lb)
            self.c1 = self.backspan.eisx(self.Lb)

        self.c2 = 0
        self.c3 = self.c1
        self.c4 = self.c1*self.a + self.c2 - self.c3*a
        self.c5 = 0.5*self.w_tot*self.b**2 + self.ml*self.b - (1.0/6.0)*self.w1*(self.b-self.a)**3 + self.c3
        self.c6 = (1.0/6.0)*self.w_tot*self.b**3 + 0.5*self.ml*self.b**2 - (1.0/24.0)*self.w1*(self.b-self.a)**4 + self.c3*self.b + self.c4 - self.c5*self.b

        arrow_height = self.w1/12.0
        #30 degree arrow
        arrow_plus_start= self.a+(arrow_height*math.tan(math.radians(30)))
        arrow_minus_start= self.a-(arrow_height*math.tan(math.radians(30)))
        arrow_plus_end= self.b+(arrow_height*math.tan(math.radians(30)))
        arrow_minus_end= self.b-(arrow_height*math.tan(math.radians(30)))

        self.x_graph=[arrow_minus_start,self.a,arrow_plus_start,self.a,self.a,self.b,self.b,arrow_minus_end,self.b,arrow_plus_end]
        self.y_graph=[arrow_height,0,arrow_height,0,self.w1,self.w1,0,arrow_height,0,arrow_height]

    def chart_load(self, x_scale=0, y_scale=0, arrows=0):
        if arrows == 1:
            arrow_height = self.w1/12.0
            #30 degree arrow
            arrow_plus_start= self.a+(arrow_height*math.tan(math.radians(30)))
            arrow_minus_start= self.a-(arrow_height*math.tan(math.radians(30)))
            arrow_plus_end= self.b+(arrow_height*math.tan(math.radians(30)))
            arrow_minus_end= self.b-(arrow_height*math.tan(math.radians(30)))

            x=[arrow_minus_start,self.a,arrow_plus_start,self.a,self.a,self.b,self.b,arrow_minus_end,self.b,arrow_plus_end]
            x = [i*x_scale for i in x]
            y=[arrow_height,0,arrow_height,0,self.w1,self.w1,0,arrow_height,0,arrow_height]
            y = [j*y_scale for j in y]
        else:
            x=[self.a,self.a,self.b,self.b]
            x = [i*x_scale for i in x]
            y=[0,self.w1,self.w1,0]
            y = [j*y_scale for j in y]

        return x,y

    def piece_functions(self):
        '''
        Returns the general piecwise function in the form of two lists
        # list1 is the polynomial coeficients of order [c0,c1x,c2x^2,...,cnx^n]
        # where the list values will only by the cn's*
        # list 2 will be the range over which the function piece applies
        # 0 <= a would be [0,a] **note it will be assumed the the eqality is <= not <
        # rerturned lists will be [[[list11],[list21]],....,[[list1n],[list2n]]
        # where n is the total number of functions to capture the range from
        # 0 to the full span, L of the beam
        '''

        v = ([
            [[self.rl],[0,self.a]],
            [[self.rl+(self.w1*self.a),-self.w1],[self.a,self.b]],
            [[0],[self.b,self.L]]
            ])

        m = ([
            [[self.ml, self.rl],[0,self.a]],
            [[self.ml-(0.5*self.a*self.a*self.w1),self.rl+(self.a*self.w1),-0.5*self.w1],[self.a,self.b]],
            [[0],[self.b,self.L]]
            ])

        eis = ([
                [[self.c1,self.ml,0.5*self.rl],[0,self.a]],
                [[self.c3+((1.0/6.0)*self.a*self.a*self.a*self.w1),
                self.ml-(0.5*self.a*self.a*self.w1),
                (0.5*self.rl)+(0.5*self.a*self.w1),
                ((-1.0/6.0)*self.w1)],[self.a,self.b]],
                [[self.c5],[self.b,self.L]]
                ])

        eid = ([
                # Range 0 to a
                [[self.c2,self.c1,0.5*self.ml,(1.0/6.0)*self.rl],[0,self.a]],
                # Range a to b
                [[self.c4-((1.0/24.0)*math.pow(self.a,4)*self.w1),      #x^0
                ((1.0/6.0)*math.pow(self.a,3)*self.w1)+self.c3,         #x^1
                ((-0.25)*math.pow(self.a,2)*self.w1)+ (0.5*self.ml),    #x^2
                ((1.0/6.0)*self.a*self.w1)+ ((1.0/6.0)*self.rl),        #x^3
                ((-1.0/24.0)*self.w1)],[self.a,self.b]],               #x^4
                # Range b to L
                [[self.c6,self.c5],[self.b,self.L]]
                ])

        vs = PieceFunctionString(v)
        ms = PieceFunctionString(m)
        eiss = PieceFunctionString(eis)
        eids = PieceFunctionString(eid)

        return [v,m,eis,eid],[vs,ms,eiss,eids]

    def fef(self):
        # Fixed End Forces
        RL = self.rl
        RR = 0
        ML = self.ml
        MR = 0

        return [RL,ML,RR,MR]

    def v(self,x):
        iters = len(x)
        v=zeros(iters)

        for i in range(0,iters):
            if x[i] <= self.a:
                v[i] = self.rl
            elif x[i]<=self.b:
                v[i] = self.rl - self.w1*(x[i]-self.a)
            else:
                v[i] = 0
        return v

    def m(self,x):
        iters = len(x)
        m=zeros(iters)

        for i in range(0,iters):
            if x[i] <= self.a:
                m[i] = self.rl*x[i] + self.ml
            elif x[i] <= self.b:
                m[i] = self.rl*x[i] + self.ml - (self.w1*(x[i]-self.a)*((x[i]-self.a)/2))
            else:
                m[i] = 0
        return m

    def eis(self,x):
        iters = len(x)
        eis=zeros(iters)

        for i in range(0,iters):
            if x[i] <= self.a:
                eis[i] = 0.5*self.rl*x[i]**2 + self.ml*x[i] + self.c1
            elif x[i] <= self.b:
                eis[i] = 0.5*self.rl*x[i]**2 + self.ml*x[i] - ((1.0/6.0) * self.w1 * (x[i]-self.a)**3) + self.c3
            else:
                eis[i] = self.c5
        return eis

    def eid(self,x):
        iters = len(x)
        eid=zeros(iters)

        for i in range(0,iters):
            if x[i] <= self.a:
                eid[i] = ((1.0/6.0)*self.rl*x[i]*x[i]*x[i]+
                            0.5*self.ml*x[i]*x[i] +
                            self.c1 * x[i] +
                            self.c2)
            elif x[i] <= self.b:
                eid[i] = ((1.0/6.0)*self.rl*x[i]*x[i]*x[i] +
                            0.5*self.ml*x[i]*x[i] -
                            ((1.0/24.0)*self.w1*(x[i]-self.a)**4) +
                            self.c3*x[i] +
                            self.c4)
            else:
                eid[i] = self.c5*x[i] + self.c6
        return eid

    def vx(self,x):
        x = float(x)
        if x <= self.a:
            v = self.w_tot
        elif x<=self.b:
            v = self.w_tot - self.w1*(x-self.a)
        else:
            v = 0
        return v

    def mx(self,x):
        x = float(x)
        if x <= self.a:
            m = self.rl*x + self.ml
        elif x <= self.b:
            m = self.rl*x + self.ml - (self.w1*(x-self.a)*((x-self.a)/2))
        else:
            m = 0
        return m

    def eisx(self,x):
        if x <= self.a:
            eis = 0.5*self.rl*x**2 + self.ml*x + self.c1
        elif x <= self.b:
            eis = 0.5*self.rl*x**2 + self.ml*x - ((1.0/6.0) * self.w1 * (x-self.a)**3) + self.c3
        else:
            eis = self.c5
        return eis

    def eidx(self,x):
        if x <= self.a:
            eid = (1.0/6.0)*self.rl*x**2 + 0.5*self.ml*x**2 + self.c1 * x + self.c2
        elif x <= self.b:
            eid = (1.0/6.0)*self.rl*x**3 + 0.5*self.ml*x**2 - (1.0/24.0)*self.w1*(x-self.a)**4 + self.c3*x + self.c4
        else:
            eid = self.c5*x + self.c6
        return eid

class cant_right_trap:
    def __init__(self, w1, w2, a, b, L, Lb):

        self.w1 = float(w1)
        self.w2 = float(w2)
        self.a = float(a)
        self.L = float(L)
        self.b = float(b)
        self.Lb = float(Lb)
        self.c = self.b-self.a

        self.kind = 'TRAP'

        self.error = ''

        if self.a > self.b:
            self.error = 'Error a > b'
            self.error = 'Error a > b'
        elif self.a > self.L:
            self.error = 'Error a > l'
            self.error = 'Error a > l'
        elif self.b > self.L:
            self.error = 'Error b > l'
            self.error = 'Error b > l'
        elif sign(self.w1) != sign(self.w2) and self.w1 !=0 and self.w2 !=0:
            self.error = 'Error w1 and w2 change direction'
            self.error = 'Error w1 and w2 change direction'
        else:
            pass

        self.w = 0.5*(self.w1+self.w2)*self.c
        self.d = self.a+(((self.w1+(2*self.w2))/(3*(self.w2+self.w1)))*self.c)
        self.s = (self.w1-self.w2)/self.c
        self.rl = self.w
        self.rr = 0
        self.ml = -1*self.w*self.d

        # 0 length backspan indicates fixed-free beam initialize slope to 0
        if Lb == 0:
            self.backspan = no_load(0)
            self.c1 = 0
        else:
            self.backspan = point_moment(-1.0*self.ml,self.Lb,self.Lb)
            self.c1 = self.backspan.eisx(self.Lb)

        self.c2 = 0

        self.c3 = self.ml - (1.0/6.0)*self.s*self.a**3 + 0.5*(self.s*self.a + self.w1)*self.a**2 - 0.5*(self.s*self.a + 2*self.w1)*self.a**2

        self.c4 = self.c1 - (1.0/24.0)*self.s*self.a**4 + (1.0/6.0)*((self.s*self.a)+self.w1)*self.a**3 - 0.25*((self.s*self.a)+(2*self.w1))*self.a**3 - self.c3*self.a + self.ml*self.a
        self.c5 = self.c1*self.a + self.c2 - self.c4*self.a - (1.0/120.0)*self.s*self.a**5 + (1.0/24.0)*((self.s*self.a)+self.w1)*self.a**4 - (1.0/12.0)*((self.s*self.a)+(2*self.w1))*self.a**4 + 0.5*self.ml*self.a**2 - 0.5*self.c3*self.a**2

        self.c6 = (0.5*self.rl*self.b**2)+self.c3*self.b + (1.0/24.0)*self.s*self.b**4 - (1.0/6.0)*((self.s*self.a)+self.w1)*self.b**3 + 0.25*((self.s*self.a)+(2*self.w1))*self.a*self.b**2 + self.c4
        self.c7 = ((1.0/6.0)*self.rl*self.b**3) + 0.5*self.c3*self.b**2 + (1.0/120.0)*self.s*self.b**5 - (1.0/24.0)*((self.s*self.a)+self.w1)*self.b**4 + (1.0/12.0)*((self.s*self.a)+(2*self.w1))*self.a*self.b**3 + self.c4*self.b + self.c5 - self.c6*self.b

        arrow_height = self.w1/6.0
        arrow_height2 = self.w2/6.0
        #30 degree arrow
        arrow_plus_start= self.a+(arrow_height*math.tan(math.radians(30)))
        arrow_minus_start= self.a-(arrow_height*math.tan(math.radians(30)))
        arrow_plus_end= self.b+(arrow_height2*math.tan(math.radians(30)))
        arrow_minus_end= self.b-(arrow_height2*math.tan(math.radians(30)))

        self.x_graph=[arrow_minus_start,self.a,arrow_plus_start,self.a,self.a,self.b,self.b,arrow_minus_end,self.b,arrow_plus_end]
        self.y_graph=[arrow_height,0,arrow_height,0,self.w1,self.w2,0,arrow_height2,0,arrow_height2]

    def chart_load(self, x_scale=0, y_scale=0, arrows=0):
        if arrows == 1:
            arrow_height = self.w1/6.0
            arrow_height2 = self.w2/6.0
            #30 degree arrow
            arrow_plus_start= self.a+(arrow_height*math.tan(math.radians(30)))
            arrow_minus_start= self.a-(arrow_height*math.tan(math.radians(30)))
            arrow_plus_end= self.b+(arrow_height2*math.tan(math.radians(30)))
            arrow_minus_end= self.b-(arrow_height2*math.tan(math.radians(30)))

            x=[arrow_minus_start,self.a,arrow_plus_start,self.a,self.a,self.b,self.b,arrow_minus_end,self.b,arrow_plus_end]
            x = [i*x_scale for i in x]
            y=[arrow_height,0,arrow_height,0,self.w1,self.w2,0,arrow_height2,0,arrow_height2]
            y = [j*y_scale for j in y]
        else:
            x=[self.a,self.a,self.b,self.b]
            x = [i*x_scale for i in x]
            y=[0,self.w1,self.w2,0]
            y = [j*y_scale for j in y]

        return x,y

    def piece_functions(self):
        '''
        Returns the general piecwise function in the form of two lists
        # list1 is the polynomial coeficients of order [c0,c1x,c2x^2,...,cnx^n]
        # where the list values will only by the cn's*
        # list 2 will be the range over which the function piece applies
        # 0 <= a would be [0,a] **note it will be assumed the the eqality is <= not <
        # rerturned lists will be [[[list11],[list21]],....,[[list1n],[list2n]]
        # where n is the total number of functions to capture the range from
        # 0 to the full span, L of the beam
        '''

        v = ([
            # Range 0 to a
            [[self.rl],[0,self.a]],
            # Range a to b
            [[(0.5*math.pow(self.a,2)*self.s) + (self.a*self.w1) + self.rl,     #x^0
            (-1.0*self.w1) - (self.a*self.s),                                   #x^1
            0.5*self.s],                                                        #x^2
            [self.a,self.b]],
            # Range b to L
            [[0],[self.b,self.L]]
            ])

        m = ([
            # Range 0 to a
            [[self.ml, self.rl],[0,self.a]],
            # Range a to b
            [[self.c3,                                                      #x^0
            (0.5*math.pow(self.a,2)*self.s)+ (self.a*self.w1) + self.rl,    #x^1
            (-0.5*self.a*self.s)-(0.5*self.w1),                             #x^2
            (1/6.0)*self.s],                                                #x^3
            [self.a,self.b]],
            # Range b to L
            [[0],[self.b,self.L]]
            ])

        eis = ([
                # Range 0 to a
                [[self.c1,self.ml,0.5*self.rl],[0,self.a]],
                # Range a to b
                [[self.c4,#x^0
                self.c3,#x^1
                (0.25*math.pow(self.a,2)*self.s)+(0.5*self.a*self.w1)+(0.5*self.rl),#x^2
                ((-1/6.0)*self.a*self.s) - ((1/6.0)*self.w1),#x^3
                (1/24.0)*self.s],#x^4
                [self.a,self.b]],
                # Range b to L
                [[self.c6],[self.b,self.L]]
                ])

        eid = ([
                # Range 0 to a
                [[self.c2,#x^0
                self.c1,#x^1
                0.5*self.ml,#x^2
                ((1.0/6.0)*self.rl),#x^3
                ],
                [0,self.a]],
                # Range a to b
                [[self.c5,#x^0
                self.c4,#x^1
                0.5*self.c3,#x^2
                ((1/12.0)*math.pow(self.a,2)*self.s)+
                ((1/6.0)*self.a*self.w1) + ((1/6.0)*self.rl),#x^3
                ((-1/24.0)*self.a*self.s) - ((1/24.0)*self.w1),#x^4
                (1/120.0)*self.s],#x^5
                [self.a,self.b]],
                # Range b to L
                [[self.c7,self.c6],[self.b,self.L]]
                ])

        vs = PieceFunctionString(v)
        ms = PieceFunctionString(m)
        eiss = PieceFunctionString(eis)
        eids = PieceFunctionString(eid)

        return [v,m,eis,eid],[vs,ms,eiss,eids]

    def fef(self):
        # Fixed End Forces
        RL = self.rl
        RR = 0
        ML = self.ml
        MR = 0

        return [RL,ML,RR,MR]

    def v(self,x):
        iters = len(x)
        v=zeros(iters)

        for i in range(0,iters):
            if x[i] <= self.a:
                v[i] = self.rl
            elif x[i]<=self.b:
                v[i] = self.rl + 0.5*self.s*x[i]**2 - x[i]*((self.s*self.a)+self.w1) + 0.5*self.a*((self.s*self.a)+(2*self.w1))
            else:
                v[i] = 0
        return v

    def m(self,x):
        iters = len(x)
        m=zeros(iters)

        for i in range(0,iters):
            if x[i] <= self.a:
                m[i] = self.rl*x[i] + self.ml
            elif x[i] <= self.b:
                m[i] = self.rl*x[i] + self.c3 + (1.0/6.0)*self.s*x[i]**3 - 0.5*((self.s*self.a)+self.w1)*x[i]**2 + 0.5*((self.s*self.a)+(2*self.w1))*self.a*x[i]
            else:
                m[i] = 0
        return m

    def eis(self,x):
        iters = len(x)
        eis=zeros(iters)

        for i in range(0,iters):
            if x[i] <= self.a:
                eis[i] = (0.5*self.rl*x[i]**2)+self.ml*x[i]+self.c1
            elif x[i] <= self.b:
                eis[i] = (0.5*self.rl*x[i]**2)+self.c3*x[i] + (1.0/24.0)*self.s*x[i]**4 - (1.0/6.0)*((self.s*self.a)+self.w1)*x[i]**3 + 0.25*((self.s*self.a)+(2*self.w1))*self.a*x[i]**2 + self.c4
            else:
                eis[i] = self.c6
        return eis

    def eid(self,x):
        iters = len(x)
        eid=zeros(iters)

        for i in range(0,iters):
            if x[i] <= self.a:
                eid[i] = ((1.0/6.0)*self.rl*x[i]**3)+ 0.5*self.ml*x[i]**2 + self.c1*x[i] + self.c2
            elif x[i] <= self.b:
                eid[i] = ((1.0/6.0)*self.rl*x[i]**3) + 0.5*self.c3*x[i]**2 + (1.0/120.0)*self.s*x[i]**5 - (1.0/24.0)*((self.s*self.a)+self.w1)*x[i]**4 + (1.0/12.0)*((self.s*self.a)+(2*self.w1))*self.a*x[i]**3 + self.c4*x[i] + self.c5
            else:
                eid[i] = self.c6*x[i] + self.c7
        return eid

    def vx(self,x):
        if x <= self.a:
            v= self.rl
        elif x<=self.b:
            v= self.rl + 0.5*self.s*x**2 - x*((self.s*self.a)+self.w1) + 0.5*self.a*((self.s*self.a)+(2*self.w1))
        else:
            v =0
        return v

    def mx(self,x):
        if x <= self.a:
            m = self.rl*x + self.ml
        elif x <= self.b:
            m = self.rl*x + self.c3 + (1.0/6.0)*self.s*x**3 - 0.5*((self.s*self.a)+self.w1)*x**2 + 0.5*((self.s*self.a)+(2*self.w1))*self.a*x
        else:
            m = 0
        return m

    def eisx(self,x):
        if x <= self.a:
            eis = (0.5*self.rl*x**2)+self.ml*x+self.c1
        elif x <= self.b:
            eis = (0.5*self.rl*x**2)+self.c3*x + (1.0/24.0)*self.s*x**4 - (1.0/6.0)*((self.s*self.a)+self.w1)*x**3 + 0.25*((self.s*self.a)+(2*self.w1))*self.a*x**2 + self.c4
        else:
            eis = self.c6
        return eis

    def eidx(self,x):
        if x <= self.a:
            eid = ((1.0/6.0)*self.rl*x**3)+ 0.5*self.ml*x**2 + self.c1*x + self.c2
        elif x <= self.b:
            eid = ((1.0/6.0)*self.rl*x**3) + 0.5*self.c3*x**2 + (1.0/120.0)*self.s*x**5 - (1.0/24.0)*((self.s*self.a)+self.w1)*x**4 + (1.0/12.0)*((self.s*self.a)+(2*self.w1))*self.a*x**3 + self.c4*x + self.c5
        else:
            eid = self.c6*x + self.c7
        return eid

class cant_left_nl:
    def __init__(self, slope, L):
        self.L = float(L)
        self.slope = float(slope)
        self.c1 = self.slope
        self.c2 = -1.0*self.c1*self.L

        self.kind = 'NL'

        self.rr = 0
        self.rl = 0
        self.mr = 0

        self.x_graph = [0]
        self.y_graph = [0]

    def chart_load(self, x_scale=0, y_scale=0, arrows=0):
        x = [0]
        y = [0]

        return x,y

    def piece_functions(self):
        '''
        Returns the general piecwise function in the form of two lists
        # list1 is the polynomial coeficients of order [c0,c1x,c2x^2,...,cnx^n]
        # where the list values will only by the cn's*
        # list 2 will be the range over which the function piece applies
        # 0 <= a would be [0,a] **note it will be assumed the the eqality is <= not <
        # rerturned lists will be [[[list11],[list21]],....,[[list1n],[list2n]]
        # where n is the total number of functions to capture the range from
        # 0 to the full span, L of the beam
        '''

        v = [[[0],[0,self.L]]]

        m = [[[0],[0,self.L]]]

        eis = [[[self.c1],[0,self.L]]]

        eid = [[[self.c2, self.c1],[0,self.L]]]

        vs = PieceFunctionString(v)
        ms = PieceFunctionString(m)
        eiss = PieceFunctionString(eis)
        eids = PieceFunctionString(eid)

        return [v,m,eis,eid],[vs,ms,eiss,eids]

    def fef(self):
        # Fixed End Forces
        RL = 0
        RR = 0
        ML = 0
        MR = 0

        return [RL,ML,RR,MR]

    def v(self,x):
        iters = len(x)
        v=zeros(iters)

        return v

    def m(self,x):
        iters = len(x)
        m=zeros(iters)

        return m

    def eis(self,x):
        iters = len(x)
        eis=zeros(iters)
        for i in range(0,iters):
            eis[i] = self.c1

        return eis

    def eid(self,x):
        iters = len(x)
        eid=zeros(iters)
        for i in range(0,iters):
            eid[i] = self.c1* x[i] + self.c2

        return eid

    def vx(self,x):
        v=0

        return v

    def mx(self,x):
        m=0

        return m

    def eisx(self,x):
        eis = self.c1

        return eis

    def eidx(self,x):
        eid = self.c1 * x + self.c2

        return eid

class cant_left_point:
    def __init__(self, p, a, L,Lb):

        self.p = float(p)
        self.a = float(a)
        self.L = float(L)
        self.Lb = float(Lb)

        self.kind = 'Point'

        self.error = ''

        if self.a > self.L:
            self.error = 'Error a > l'
            self.error = 'Error a > l'

        self.rr = self.p
        self.rl = 0
        self.mr = -1*self.p*(self.L-self.a)

        # 0 length backspan indicates fixed-free beam initialize slope to 0
        if self.Lb == 0:
            self.backspan = no_load(0)
            self.c3 = 0 + (0.5*self.p * (self.L-self.a)**2)
        else:
            self.backspan = point_moment(self.mr,0,self.Lb)
            self.c3 = self.backspan.eisx(0) + (0.5*self.p * (self.L-self.a)**2)

        self.c4 = ((1/6.0)*self.p*(self.L-self.a)**3) - (self.c3*self.L)
        self.c1 = self.c3
        self.c2 = (self.c3*self.a) + self.c4 - (self.c1*self.a)

        arrow_height = self.p/6.0
        #30 degree arrow
        arrow_plus= self.a+(arrow_height*math.tan(math.radians(30)))
        arrow_minus= self.a-(arrow_height*math.tan(math.radians(30)))

        self.x_graph=[arrow_minus,self.a,arrow_plus,self.a,self.a]
        self.y_graph=[arrow_height,0,arrow_height,0,self.p]

    def chart_load(self, x_scale=0, y_scale=0, arrows=0):
        if arrows == 1:
            arrow_height = self.p/6.0
            #30 degree arrow
            arrow_plus= self.a+(arrow_height*math.tan(math.radians(30)))
            arrow_minus= self.a-(arrow_height*math.tan(math.radians(30)))

            x=[arrow_minus,self.a,arrow_plus,self.a,self.a]
            x = [i*x_scale for i in x]
            y=[arrow_height,0,arrow_height,0,self.p]
            y = [j*y_scale for j in y]
        else:
            x=[self.a,self.a]
            x = [i*x_scale for i in x]
            y=[0,self.p]
            y = [j*y_scale for j in y]

        return x,y

    def piece_functions(self):
        '''
        Returns the general piecwise function in the form of two lists
        # list1 is the polynomial coeficients of order [c0,c1x,c2x^2,...,cnx^n]
        # where the list values will only by the cn's*
        # list 2 will be the range over which the function piece applies
        # 0 <= a would be [0,a] **note it will be assumed the the eqality is <= not <
        # rerturned lists will be [[[list11],[list21]],....,[[list1n],[list2n]]
        # where n is the total number of functions to capture the range from
        # 0 to the full span, L of the beam
        '''


        v = [[[0],[0,self.a]],[[-1.0*self.p],[self.a,self.L]]]

        m = [[[0],[0,self.a]],[[self.p*self.a,-1.0*self.p],[self.a,self.L]]]

        eis = [[[self.c1],[0,self.a]],[[-0.5*self.a*self.a*self.p+self.c3,self.a*self.p, -0.5*self.p],[self.a,self.L]]]

        eid = [[[self.c2,self.c1],[0,self.a]],[[self.c4+((self.a*self.a*self.a*self.p)*(1/6.0)), self.c3-(0.5*self.a*self.a*self.p),0.5*self.a*self.p,(-1/6.0)*self.p],[self.a,self.L]]]

        vs = PieceFunctionString(v)
        ms = PieceFunctionString(m)
        eiss = PieceFunctionString(eis)
        eids = PieceFunctionString(eid)

        return [v,m,eis,eid],[vs,ms,eiss,eids]

    def fef(self):
        # Fixed End Forces
        RL = 0
        RR = self.rr
        ML = 0
        MR = self.mr

        return [RL,ML,RR,MR]

    def v(self,x):
        iters = len(x)
        v=zeros(iters)

        for i in range(0,iters):
            if x[i]<=self.a:
                v[i] = 0
            else:
                v[i] = -1*self.p
        return v

    def m(self,x):
        iters = len(x)
        m=zeros(iters)

        for i in range(0,iters):
            if x[i]<=self.a:
                m[i] = 0
            else:
                m[i] = -1*self.p * (x[i] - self.a)
        return m

    def eis(self,x):
        iters = len(x)
        eis=zeros(iters)

        for i in range(0,iters):
            if x[i]<=self.a:
                eis[i] = self.c1
            else:
                eis[i] = (-0.5*self.p * (x[i]-self.a)**2) + self.c3
        return eis

    def eid(self, x):
        iters = len(x)
        eid=zeros(iters)

        for i in range(0,iters):
            if x[i]<=self.a:
                eid[i] = self.c1*x[i] + self.c2
            else:
                eid[i] = (-1/6.0)*self.p*(x[i]-self.a)**3 + self.c3*x[i] + self.c4
        return eid

    def vx(self,x):
        if x<=self.a:
            v = 0
        else:
            v = -1*self.p
        return v

    def mx(self,x):
        if x<=self.a:
            m = 0
        else:
            m = -1*self.p * (x - self.a)
        return m

    def eisx(self,x):
        if x<=self.a:
            eis = self.c1
        else:
            eis  = (-0.5*self.p * (x-self.a)**2) + self.c3
        return eis

    def eidx(self, x):
        if x<=self.a:
            eid = self.c1*x + self.c2
        else:
            eid = (-1/6.0)*self.p*(x-self.a)**3 + self.c3*x + self.c4

        return eid

class cant_left_point_moment:
    def __init__(self, ma, a, L,Lb):

        self.ma = float(ma)
        self.a = float(a)
        self.L = float(L)
        self.Lb = float(Lb)

        self.kind = 'Moment'

        self.error = ''

        if self.a > self.L:
            self.error = 'Error a > l'
            self.error = 'Error a > l'

        self.rr = 0
        self.rl = 0
        self.mr = self.ma

        # 0 length backspan indicates fixed-free beam initialize slope to 0
        if Lb == 0:
            self.backspan = no_load(0)
            self.c3 = 0 - (self.ma*self.L)
        else:
            self.backspan = point_moment(self.mr,0,Lb)
            self.c3 = self.backspan.eisx(0) - (self.ma*self.L)

        self.c4 = (-0.5*self.ma*self.L**2) - self.c3*self.L
        self.c1 = (1.0*self.ma*self.a) + self.c3
        self.c2 = 0.5*self.ma*self.a**2 + self.c3*self.a + self.c4 - self.c1*self.a

        r = (self.ma/2.0)
        arrow_height = r/6.0
        #30 degree arrow
        arrow_minus= (arrow_height*math.tan(math.radians(30)))

        if self.ma <0:
            self.x_graph = [self.a,self.a,self.a]
            self.y_graph = [r,0,-r]
            x=0
            y=0
            for a in range(-90, 181):
                x = self.a+(r*math.cos(math.radians(a)))
                y = 0+(r*math.sin(math.radians(a)))
                self.x_graph.append(x)
                self.y_graph.append(y)

            self.x_graph.append(x-arrow_minus)
            self.y_graph.append(y+arrow_height)
            self.x_graph.append(x)
            self.y_graph.append(y)
            self.x_graph.append(x+arrow_minus)
            self.y_graph.append(y+arrow_height)
        else:
            self.x_graph = [self.a-r,self.a,self.a+r, self.a+r-arrow_minus,self.a+r,self.a+r+arrow_minus,self.a+r]
            self.y_graph = [0,0,0,arrow_height,0,arrow_height,0]
            x=0
            y=0
            for a in range(0,271):
                x = self.a+(r*math.cos(math.radians(a)))
                y = 0+(r*math.sin(math.radians(a)))
                self.x_graph.append(x)
                self.y_graph.append(y)

    def chart_load(self, x_scale=0, y_scale=0, arrows=0):
        r = (self.ma/2.0)

        if arrows == 1:
            arrow_height = r/6.0
            #30 degree arrow
            arrow_minus= (arrow_height*math.tan(math.radians(30)))

            if self.ma <0:
                x= [self.a,self.a,self.a]
                y = [r,0,-r]
                xi=0
                yi=0
                for a in range(-90, 181):
                    xi = self.a+(r*math.cos(math.radians(a)))
                    yi = 0+(r*math.sin(math.radians(a)))
                    x.append(xi)
                    y.append(yi)

                x.append(xi-arrow_minus)
                y.append(yi+arrow_height)
                x.append(xi)
                y.append(yi)
                x.append(xi+arrow_minus)
                y.append(yi+arrow_height)
            else:
                x= [self.a-r,self.a,self.a+r, self.a+r-arrow_minus,self.a+r,self.a+r+arrow_minus,self.a+r]
                y = [0,0,0,arrow_height,0,arrow_height,0]
                xi=0
                yi=0
                for a in range(0,271):
                    xi = self.a+(r*math.cos(math.radians(a)))
                    yi = 0+(r*math.sin(math.radians(a)))
                    x.append(xi)
                    y.append(yi)

        x = [i*x_scale for i in x]
        y = [j*y_scale for j in y]

        return x,y

    def piece_functions(self):
        '''
        Returns the general piecwise function in the form of two lists
        # list1 is the polynomial coeficients of order [c0,c1x,c2x^2,...,cnx^n]
        # where the list values will only by the cn's*
        # list 2 will be the range over which the function piece applies
        # 0 <= a would be [0,a] **note it will be assumed the the eqality is <= not <
        # rerturned lists will be [[[list11],[list21]],....,[[list1n],[list2n]]
        # where n is the total number of functions to capture the range from
        # 0 to the full span, L of the beam
        '''

        v = [[[0],[0,self.L]]]

        m = [[[0],[0,self.a]],[[self.ma],[self.a,self.L]]]

        eis = [[[self.c1],[0,self.a]],[[self.c3,self.ma],[self.a,self.L]]]

        eid = [[[self.c2, self.c1],[0,self.a]],[[self.c4,self.c3,0.5*self.ma],[self.a,self.L]]]

        vs = PieceFunctionString(v)
        ms = PieceFunctionString(m)
        eiss = PieceFunctionString(eis)
        eids = PieceFunctionString(eid)

        return [v,m,eis,eid],[vs,ms,eiss,eids]

    def fef(self):
        # Fixed End Forces
        RL = 0
        RR = self.rr
        ML = 0
        MR = self.mr

        return [RL,ML,RR,MR]

    def v(self,x):
        iters = len(x)
        v=zeros(iters)

        for i in range(0,iters):
            if x[i]<=self.a:
                v[i] = 0
            else:
                v[i] = 0
        return v

    def m(self,x):
        iters = len(x)
        m=zeros(iters)

        for i in range(0,iters):
            if x[i]<=self.a:
                m[i] = 0
            else:
                m[i] = self.ma
        return m

    def eis(self,x):
        iters = len(x)
        eis=zeros(iters)

        for i in range(0,iters):
            if x[i]<=self.a:
                eis[i] = self.c1
            else:
                eis[i] = (self.ma * x[i]) + self.c3
        return eis

    def eid(self, x):
        iters = len(x)
        eid=zeros(iters)

        for i in range(0,iters):
            if x[i]<=self.a:
                eid[i] = self.c1*x[i] + self.c2
            else:
                eid[i] = (0.5)*self.ma*x[i]**2 + self.c3*x[i] + self.c4
        return eid

    def vx(self,x):
        if x<=self.a:
            v = 0
        else:
            v = 0
        return v

    def mx(self,x):
        if x<=self.a:
            m = 0
        else:
            m = self.ma
        return m

    def eisx(self,x):
        if x<=self.a:
            eis = self.c1
        else:
            eis = (self.ma * x) + self.c3
        return eis

    def eidx(self, x):
        if x<=self.a:
            eid = self.c1*x + self.c2
        else:
            eid = (0.5)*self.ma*x**2 + self.c3*x + self.c4
        return eid

class cant_left_udl:
    def __init__(self, w1, a, b, L, Lb):

        self.w1 = float(w1)
        self.a = float(a)
        self.L = float(L)
        self.Lb = float(Lb)
        self.b = float(b)
        self.c = self.b-self.a
        self.w_tot = self.w1*self.c

        self.kind = 'UDL'

        self.error = ''

        if self.a > self.b:
            self.error = 'Error a > b'
            self.error = 'Error a > b'
        elif self.a > self.L:
            self.error = 'Error a > l'
            self.error = 'Error a > l'
        elif self.b > self.L:
            self.error = 'Error b > l'
            self.error = 'Error b > l'
        else:
            pass

        self.rr = self.w_tot
        self.rl = 0
        self.mr = -1.0*self.w_tot*(self.L-(a+(self.c/2.0)))

        # 0 length backspan indicates fixed-free beam initialize slope to 0
        if Lb == 0:
            self.backspan = no_load(0)
            self.c5 = 0 + (0.5 * self.w_tot * (self.L - (self.a + (0.5*self.c)))**2)
        else:
            self.backspan = point_moment(self.mr,0,Lb)
            self.c5 = self.backspan.eisx(0) + (0.5 * self.w_tot * (self.L - (self.a + (0.5*self.c)))**2)

        self.c6 = ((1.0/6.0)*self.w_tot * (self.L - (self.a + (0.5*self.c)))**3) - (self.c5*self.L)
        self.c3 =((-0.5)*self.w_tot * (self.b - (self.a + (0.5*self.c)))**2) + self.c5 + ((1.0/6.0)*self.w1*(b-a)**3)
        self.c1 = self.c3
        self.c4 = ((-1.0/6.0)*self.w_tot * (self.b - (self.a + (0.5*self.c)))**3) + (self.c5*self.b) + self.c6 + ((1.0/24.0)*self.w1*(self.b-self.a)**4) - (self.c3*self.b)
        self.c2 = (self.c3*self.a) + self.c4 - (self.c1*self.a)

        arrow_height = self.w1/12.0
        #30 degree arrow
        arrow_plus_start= self.a+(arrow_height*math.tan(math.radians(30)))
        arrow_minus_start= self.a-(arrow_height*math.tan(math.radians(30)))
        arrow_plus_end= self.b+(arrow_height*math.tan(math.radians(30)))
        arrow_minus_end= self.b-(arrow_height*math.tan(math.radians(30)))

        self.x_graph=[arrow_minus_start,self.a,arrow_plus_start,self.a,self.a,self.b,self.b,arrow_minus_end,self.b,arrow_plus_end]
        self.y_graph=[arrow_height,0,arrow_height,0,self.w1,self.w1,0,arrow_height,0,arrow_height]

    def chart_load(self, x_scale=0, y_scale=0, arrows=0):
        if arrows == 1:
            arrow_height = self.w1/12.0
            #30 degree arrow
            arrow_plus_start= self.a+(arrow_height*math.tan(math.radians(30)))
            arrow_minus_start= self.a-(arrow_height*math.tan(math.radians(30)))
            arrow_plus_end= self.b+(arrow_height*math.tan(math.radians(30)))
            arrow_minus_end= self.b-(arrow_height*math.tan(math.radians(30)))

            x=[arrow_minus_start,self.a,arrow_plus_start,self.a,self.a,self.b,self.b,arrow_minus_end,self.b,arrow_plus_end]
            x = [i*x_scale for i in x]
            y=[arrow_height,0,arrow_height,0,self.w1,self.w1,0,arrow_height,0,arrow_height]
            y = [j*y_scale for j in y]
        else:
            x=[self.a,self.a,self.b,self.b]
            x = [i*x_scale for i in x]
            y=[0,self.w1,self.w1,0]
            y = [j*y_scale for j in y]

        return x,y

    def piece_functions(self):
        '''
        Returns the general piecwise function in the form of two lists
        # list1 is the polynomial coeficients of order [c0,c1x,c2x^2,...,cnx^n]
        # where the list values will only by the cn's*
        # list 2 will be the range over which the function piece applies
        # 0 <= a would be [0,a] **note it will be assumed the the eqality is <= not <
        # rerturned lists will be [[[list11],[list21]],....,[[list1n],[list2n]]
        # where n is the total number of functions to capture the range from
        # 0 to the full span, L of the beam
        '''

        v = ([
            [[0],[0,self.a]],
            [[self.w1*self.a,
            -1.0*self.w1],
            [self.a,self.b]],
            [[-1.0*self.w_tot],[self.b,self.L]]
            ])

        m = ([
            # Range 0 to a
            [[0],[0,self.a]],
            # Range a to b
            [[-0.5*math.pow(self.a,2)*self.w1,
            self.a*self.w1,
            -0.5*self.w1],
            [self.a,self.b]],
            # Range b to L
            [[self.a*self.w_tot + 0.5*self.c*self.w_tot,
             -1.0*self.w_tot],
             [self.b,self.L]]
            ])

        eis = ([
                # Range 0 to a
                [[self.c1],[0,self.a]],
                # Range a to b
                [[(1/6.0)*math.pow(self.a,3)*self.w1 + self.c3,#x^0
                -0.5*math.pow(self.a,2)*self.w1,#x^1
                0.5*self.a*self.w1,#x^2
                (-1/6.0)*self.w1],#x^3
                [self.a,self.b]],
                # Range b to L
                [[self.c5-(0.5*math.pow(self.a,2)*self.w_tot)-
                (0.5*self.a*self.c*self.w_tot) - ((1/8.0)*math.pow(self.c,2)*self.w_tot),#x^0
                (self.a*self.w_tot)+(0.5*self.c*self.w_tot),#x^1
                -0.5*self.w_tot],#x^2
                [self.b,self.L]]
                ])

        eid = ([
                # Range 0 to a
                [[self.c2,self.c1],[0,self.a]],
                # Range a to b
                [[self.c4-((1/24.0)*math.pow(self.a,4)*self.w1),#x^0
                (1/6.0)*math.pow(self.a,3)*self.w1+self.c3,#x^1
                -0.25*math.pow(self.a,2)*self.w1,#x^2
                (1/6.0)*self.a*self.w1,#x^3
                (-1/24.0)*self.w1],#x^4
                [self.a,self.b]],
                # Range b to L
                [[((1/6.0)*math.pow(self.a,3)*self.w_tot)+
                (0.25*math.pow(self.a,2)*self.c*self.w_tot)+
                (0.125*self.a*math.pow(self.c,2)*self.w_tot)+
                ((1/48.0)*math.pow(self.c,3)*self.w_tot)+self.c6,#x^0
                (-0.5*math.pow(self.a,2)*self.w_tot)-
                (0.5*self.a*self.c*self.w_tot)-
                (0.125*math.pow(self.c,2)*self.w_tot)+self.c5,#x^1
                (0.5*self.a*self.w_tot) + (0.25*self.c*self.w_tot),#x^2
                (-1/6.0)*self.w_tot],#x^3
                [self.b,self.L]]
                ])

        vs = PieceFunctionString(v)
        ms = PieceFunctionString(m)
        eiss = PieceFunctionString(eis)
        eids = PieceFunctionString(eid)

        return [v,m,eis,eid],[vs,ms,eiss,eids]

    def fef(self):
        # Fixed End Forces
        RL = 0
        RR = self.rr
        ML = 0
        MR = self.mr

        return [RL,ML,RR,MR]

    def v(self,x):
        iters = len(x)
        v=zeros(iters)

        for i in range(0,iters):
            if x[i] <= self.a:
                v[i] = 0
            elif x[i]<=self.b:
                v[i] = -1*self.w1*(x[i]-self.a)
            else:
                v[i] = -1*self.w_tot
        return v

    def m(self,x):
        iters = len(x)
        m=zeros(iters)

        for i in range(0,iters):
            if x[i] <= self.a:
                m[i] = 0
            elif x[i] <= self.b:
                m[i] = -0.5*self.w1*(x[i]-self.a)**2
            else:
                m[i] = -1.0 * self.w_tot * (x[i]-(self.a+(0.5*self.c)))
        return m

    def eis(self,x):
        iters = len(x)
        eis=zeros(iters)

        for i in range(0,iters):
            if x[i] <= self.a:
                eis[i] = self.c1
            elif x[i] <= self.b:
                eis[i] = (-1.0/6.0)*self.w1*(x[i]-self.a)**3 + self.c3
            else:
                eis[i] = (-0.5 * self.w_tot * (x[i]-(self.a+(0.5*self.c)))**2) + self.c5
        return eis

    def eid(self,x):
        iters = len(x)
        eid=zeros(iters)

        for i in range(0,iters):
            if x[i] <= self.a:
                eid[i] = self.c1*x[i] + self.c2
            elif x[i] <= self.b:
                eid[i] = (-1.0/24.0)*self.w1*(x[i]-self.a)**4 + self.c3*x[i] + self.c4
            else:
                eid[i] = ((-1.0/6.0) * self.w_tot * (x[i]-(self.a+(0.5*self.c)))**3) + self.c5*x[i] + self.c6
        return eid

    def vx(self,x):
        if x <= self.a:
            v = 0
        elif x<=self.b:
            v = -1*self.w1*(x-self.a)
        else:
            v = -1*self.w_tot
        return v

    def mx(self,x):
        if x <= self.a:
            m = 0
        elif x <= self.b:
            m = -0.5*self.w1*(x-self.a)**2
        else:
            m = -1.0 * self.w_tot * (x-(self.a+(0.5*self.c)))
        return m

    def eisx(self,x):
        if x <= self.a:
            eis = self.c1
        elif x <= self.b:
            eis = (-1.0/6.0)*self.w1*(x-self.a)**3 + self.c3
        else:
            eis = (-0.5 * self.w_tot * (x-(self.a+(0.5*self.c)))**2) + self.c5
        return eis

    def eidx(self,x):
        if x <= self.a:
            eid = self.c1*x+ self.c2
        elif x <= self.b:
            eid = (-1.0/24.0)*self.w1*(x-self.a)**4 + self.c3*x + self.c4
        else:
            eid = ((-1.0/6.0) * self.w_tot * (x-(self.a+(0.5*self.c)))**3) + self.c5*x + self.c6
        return eid

class cant_left_trap:

    def __init__(self, w1, w2, a, b, L, Lb):

        self.w1 = float(w1)
        self.w2 = float(w2)
        self.a = float(a)
        self.L = float(L)
        self.b = float(b)
        self.Lb = float(Lb)
        self.c = self.b-self.a

        self.kind = 'TRAP'

        self.error = ''

        if self.a > self.b:
            self.error = 'Error a > b'
            self.error = 'Error a > b'
        elif self.a > self.L:
            self.error = 'Error a > l'
            self.error = 'Error a > l'
        elif self.b > self.L:
            self.error = 'Error b > l'
            self.error = 'Error b > l'
        elif sign(self.w1) != sign(self.w2) and self.w1 !=0 and self.w2 !=0:
            self.error = 'Error w1 and w2 change direction'
            self.error = 'Error w1 and w2 change direction'
        else:
            pass

        self.w = 0.5*(self.w1+self.w2)*self.c
        self.dl = self.a+(((self.w1+(2*self.w2))/(3*(self.w2+self.w1)))*self.c)
        self.dr = self.L-self.dl
        self.s = (self.w1-self.w2)/self.c
        self.cc = (((self.w1+(2*self.w2))/(3*(self.w2+self.w1)))*self.c) + self.a
        self.rr = self.w
        self.rl=0
        self.mr = -1*self.rr*(self.L-self.cc)

        # 0 length backspan indicates fixed-free beam initialize slope to 0
        if Lb == 0:
            self.backspan = no_load(0)
            self.c6 = 0 + (0.5*self.w*(self.L-self.cc)**2)
        else:
            self.backspan = point_moment(self.mr,0,Lb)
            self.c6 = self.backspan.eisx(0) + (0.5*self.w*(self.L-self.cc)**2)

        self.c7 = ((1.0/6.0)*self.w*(self.L-self.cc)**3) - (self.c6*self.L)
        self.c3 = -1.0*((1.0/6.0)*self.a*((self.a**2 * self.s) - (3*self.a*((self.a*self.s) + self.w1)) + (3*self.a*((self.a*self.s) + (2*self.w1)))))
        self.c4 = (-0.5*self.w*(self.b-self.cc)**2) + self.c6 - (self.c3*self.b) - ((1.0/24.0)*self.b**2 *((self.b**2 * self.s) - (4*self.b*((self.a*self.s) + self.w1)) + (6*self.a*((self.a*self.s) + (2*self.w1)))))
        self.c5 = ((-1.0/6.0)*self.w*(self.b-self.cc)**3) + (self.c6*self.b)+self.c7-(0.5*self.c3*self.b**2)-(self.c4*self.b)-((1.0/120.0)*self.b**3 *((self.b**2 * self.s) - (5*self.b*((self.a*self.s) + self.w1)) + (10*self.a*((self.a*self.s) + (2*self.w1)))))
        self.c1 = ((1.0/24.0)*self.a**2 *((self.a**2 * self.s) - (4*self.a*((self.a*self.s) + self.w1)) + (6*self.a*((self.a*self.s) + (2*self.w1))))) + (self.c3*self.a) + self.c4
        self.c2 = ((1.0/120.0)*self.a**3 *((self.a**2 * self.s) - (5*self.a*((self.a*self.s) + self.w1)) + (10*self.a*((self.a*self.s) + (2*self.w1))))) + (0.5*self.c3*self.a**2) + (self.c4*self.a) + self.c5 - (self.c1*self.a)

        arrow_height = self.w1/6.0
        arrow_height2 = self.w2/6.0
        #30 degree arrow
        arrow_plus_start= self.a+(arrow_height*math.tan(math.radians(30)))
        arrow_minus_start= self.a-(arrow_height*math.tan(math.radians(30)))
        arrow_plus_end= self.b+(arrow_height2*math.tan(math.radians(30)))
        arrow_minus_end= self.b-(arrow_height2*math.tan(math.radians(30)))

        self.x_graph=[arrow_minus_start,self.a,arrow_plus_start,self.a,self.a,self.b,self.b,arrow_minus_end,self.b,arrow_plus_end]
        self.y_graph=[arrow_height,0,arrow_height,0,self.w1,self.w2,0,arrow_height2,0,arrow_height2]

    def chart_load(self, x_scale=0, y_scale=0, arrows=0):
        if arrows == 1:
            arrow_height = self.w1/6.0
            arrow_height2 = self.w2/6.0
            #30 degree arrow
            arrow_plus_start= self.a+(arrow_height*math.tan(math.radians(30)))
            arrow_minus_start= self.a-(arrow_height*math.tan(math.radians(30)))
            arrow_plus_end= self.b+(arrow_height2*math.tan(math.radians(30)))
            arrow_minus_end= self.b-(arrow_height2*math.tan(math.radians(30)))

            x=[arrow_minus_start,self.a,arrow_plus_start,self.a,self.a,self.b,self.b,arrow_minus_end,self.b,arrow_plus_end]
            x = [i*x_scale for i in x]
            y=[arrow_height,0,arrow_height,0,self.w1,self.w2,0,arrow_height2,0,arrow_height2]
            y = [j*y_scale for j in y]
        else:
            x=[self.a,self.a,self.b,self.b]
            x = [i*x_scale for i in x]
            y=[0,self.w1,self.w2,0]
            y = [j*y_scale for j in y]

        return x,y

    def piece_functions(self):
        '''
        Returns the general piecwise function in the form of two lists
        # list1 is the polynomial coeficients of order [c0,c1x,c2x^2,...,cnx^n]
        # where the list values will only by the cn's*
        # list 2 will be the range over which the function piece applies
        # 0 <= a would be [0,a] **note it will be assumed the the eqality is <= not <
        # rerturned lists will be [[[list11],[list21]],....,[[list1n],[list2n]]
        # where n is the total number of functions to capture the range from
        # 0 to the full span, L of the beam
        '''

        v = ([
            [[0],
            [0,self.a]],
            [[(0.5*math.pow(self.a,2)*self.s)+(self.a*self.w1), #x^0
            (-1.0*self.a*self.s) - self.w1,                     #x^1
            0.5*self.s],                                        #x^2
            [self.a,self.b]],
            [[-1.0*self.rr],
            [self.b,self.L]]
            ])

        m = ([
            # Range 0 to a
            [[0],
            [0,self.a]],
            # Range a to b
            [[self.c3,                                          #x^0
            (0.5*math.pow(self.a,2)*self.s)+(self.a*self.w1),   #x^1
            (-0.5*self.a*self.s) - (0.5*self.w1),               #x^2
            (1/6.0)*self.s],                                    #x^3
            [self.a,self.b]],
            # Range b to L
            [[self.w*self.cc,   #x^0
             -1.0*self.w],      #x^1
             [self.b,self.L]]
            ])

        eis = ([
                # Range 0 to a
                [[self.c1],
                [0,self.a]],
                # Range a to b
                [[self.c4,#x^0
                self.c3,#x^1
                (0.25*math.pow(self.a,2)*self.s)+(0.5*self.a*self.w1),#x^2
                ((-1/6.0)*self.a*self.s)-((1/6.0)*self.w1),#x^3
                (1/24.0)*self.s],#x^4
                [self.a,self.b]],
                # Range b to L
                [[self.c6-(0.5*math.pow(self.cc,2)*self.w),#x^0
                self.cc*self.w,#x^1
                -0.5*self.w],#x^2
                [self.b,self.L]]
                ])

        eid = ([
                # Range 0 to a
                [[self.c2,self.c1],
                [0,self.a]],
                # Range a to b
                [[self.c5,#x^0
                self.c4,#x^1
                0.5*self.c3,#x^2
                ((1/12.0)*math.pow(self.a,2)*self.s)+((1/6.0)*self.a*self.w1),#x^3
                ((-1/24.0)*self.a*self.s)-((1/24.0)*self.w1),#x^4
                (1/120.0)*self.s],#x^5

                [self.a,self.b]],
                # Range b to L
                [[self.c7+((1/6.0)*math.pow(self.cc,3)*self.w),#x^0
                self.c6-(0.5*math.pow(self.cc,2)*self.w),#x^1
                0.5*self.cc*self.w,#x^2
                (-1/6.0)*self.w],#x^3
                [self.b,self.L]]
                ])

        vs = PieceFunctionString(v)
        ms = PieceFunctionString(m)
        eiss = PieceFunctionString(eis)
        eids = PieceFunctionString(eid)

        return [v,m,eis,eid],[vs,ms,eiss,eids]

    def fef(self):
        # Fixed End Forces
        RL = 0
        RR = self.rr
        ML = 0
        MR = self.mr

        return [RL,ML,RR,MR]

    def v(self,x):
        iters = len(x)
        v=zeros(iters)

        for i in range(0,iters):
            if x[i] <= self.a:
                v[i] = 0
            elif x[i]<=self.b:
                v[i] = (-0.5*((2*self.w1)-(self.s*(x[i]-self.a))))*(x[i]-self.a)
            else:
                v[i] = -1*self.rr
        return v

    def m(self,x):
        iters = len(x)
        m=zeros(iters)

        for i in range(0,iters):
            if x[i] <= self.a:
                m[i] = 0
            elif x[i] <= self.b:
                m[i] = ((1.0/6.0)*x[i]*((x[i]**2 * self.s) - (3*x[i]*((self.a*self.s) + self.w1)) + (3*self.a*((self.a*self.s) + (2*self.w1))))) + self.c3
            else:
                m[i] = -1*self.w*(x[i]-self.cc)
        return m

    def eis(self,x):
        iters = len(x)
        eis=zeros(iters)

        for i in range(0,iters):
            if x[i] <= self.a:
                eis[i] = self.c1
            elif x[i] <= self.b:
                eis[i] = ((1.0/24.0)*x[i]**2 *((x[i]**2 * self.s) - (4*x[i]*((self.a*self.s) + self.w1)) + (6*self.a*((self.a*self.s) + (2*self.w1))))) + (self.c3 * x[i]) + self.c4
            else:
                eis[i] = (-0.5*self.w*(x[i]-self.cc)**2) + self.c6
        return eis

    def eid(self,x):
        iters = len(x)
        eid=zeros(iters)

        for i in range(0,iters):
            if x[i] <= self.a:
                eid[i] = self.c1*x[i] + self.c2
            elif x[i] <= self.b:
                eid[i] = ((1.0/120.0)*x[i]**3 *((x[i]**2 * self.s) - (5*x[i]*((self.a*self.s) + self.w1)) + (10*self.a*((self.a*self.s) + (2*self.w1))))) + (0.5*self.c3 * x[i]**2) + (self.c4*x[i]) + self.c5
            else:
                eid[i] = ((-1.0/6.0)*self.w*(x[i]-self.cc)**3) + (self.c6*x[i]) + self.c7
        return eid

    def vx(self,x):
        if x <= self.a:
            v = 0
        elif x<=self.b:
            v= (-0.5*((2*self.w1)-(self.s*(x-self.a))))*(x-self.a)
        else:
            v = -1*self.rr
        return v

    def mx(self,x):
        if x <= self.a:
            m = 0
        elif x <= self.b:
            m = ((1.0/6.0)*x*((x**2 * self.s) - (3*x*((self.a*self.s) + self.w1)) + (3*self.a*((self.a*self.s) + (2*self.w1))))) + self.c3
        else:
            m = -1*self.w*(x-self.cc)
        return m

    def eisx(self,x):
        if x <= self.a:
            eis = self.c1
        elif x <= self.b:
            eis = ((1.0/24.0)*x**2 *((x**2 * self.s) - (4*x*((self.a*self.s) + self.w1)) + (6*self.a*((self.a*self.s) + (2*self.w1))))) + (self.c3 * x) + self.c4
        else:
            eis = (-0.5*self.w*(x-self.cc)**2) + self.c6
        return eis

    def eidx(self,x):
        if x <= self.a:
            eid = self.c1*x + self.c2
        elif x <= self.b:
            eid = ((1.0/120.0)*x**3 *((x**2 * self.s) - (5*x*((self.a*self.s) + self.w1)) + (10*self.a*((self.a*self.s) + (2*self.w1))))) + (0.5*self.c3 * x**2) + (self.c4*x) + self.c5
        else:
            eid = ((-1.0/6.0)*self.w*(x-self.cc)**3) + (self.c6*x) + self.c7
        return eid
