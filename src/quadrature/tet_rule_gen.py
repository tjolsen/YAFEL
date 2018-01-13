#!/usr/bin/python3
import itertools as it


def generate_points(x_tup, weight):
    x = []
    xarr = [i for i in x_tup]
    xarr.sort()
    x_tup = tuple(xarr)
    for p in it.permutations(x_tup, len(x_tup)-1):
        x.append(p)

    x = [xi for xi in set(x)]
    return x,[weight for xi in x]

#==================================================


class QuadratureRule:

    def __init__(self,p,x,w):
        self.p = p
        self.x = x
        self.w = w

    def __str__(self):
        # C++ code generation for quadrature rule
        s = "case %d:\n" % self.p

        s += "nodes = {\n"
        for xi in self.x:
            s += '{' + ','.join([str(xi_i) for xi_i in xi]) + '},\n'

        s += '};\n'
        s += 'weights = {' + ",\n".join([str(wi) for wi in self.w]) + '};\n'
        s += 'return;'
        
        return s

#=================================================    


def gen_quadrature_rules(qr_list):
    QRs = []
    
    for qri in qr_list:
        p = qri[0]
        x = []
        w = []
        for xwi in qri[1]:
            xi,wi = generate_points(xwi[0], xwi[1])
            x.extend(xi)
            w.extend(wi)
        QRs.append(QuadratureRule(p,x,w))

    return QRs

#=================================================



# list of tuples: (polyOrder, [((a,b,c), weight)]) where (a,b,c) are
# values to be permuted to generate quadrature points
tri_quad_rules = [
    (1, [((1.0/3.0, 1.0/3.0, 1.0/3.0), 0.5)]),
    (2, [((2.0/3.0, 1.0/6.0, 1.0/6.0), 1.0/6.0)]),
    (3, [((1.0/3.0, 1.0/3.0, 1.0/3.0), -2.81250000000000000e-01),
         ((0.6, 0.2, 0.2), 2.60416666666666667E-01)])
]



tet_quad_rules = [
    (1, [( (0.25,0.25,0.25,0.25), 1.0/6.0)]),
    (2, [((0.585410196624969, 0.138196601125011, 0.138196601125011, 0.138196601125011), 0.25/6.0)]),
    (3, [
        ( (0.25,0.25,0.25,0.25), -0.8/6.0),
        ( (0.5, 1.0/6.0, 1.0/6.0, 1.0/6.0), 0.45/6.0)
    ]),
    (4, [
        ( (0.25,0.25,0.25,0.25),
          -0.013155555555555555556),
        ( (0.785714285714286, 0.071428571428571, 0.071428571428571, 0.071428571428571),
          0.007622222222222222),
        ( (0.399403576166799, 0.399403576166799, 0.100596423833201, 0.100596423833201),
          0.024888888888888889)
    ]),
    (5, [
        ( (0.25,0.25,0.25,0.25),
          0.030283678097089),
        ( (0.0, 1.0/3.0, 1.0/3.0, 1.0/3.0),
          0.006026785714286),
        ( (8.0/11.0, 1.0/11.0, 1.0/11.0, 1.0/11.0),
          0.011645249086029),
        ( (0.066550153573664, 0.066550153573664, 0.433449846426336, 0.433449846426336),
          0.010949141561386)
    ])
]

tri_QRs = gen_quadrature_rules(tri_quad_rules)
tet_QRs = gen_quadrature_rules(tet_quad_rules)
for qr in tet_QRs:
    print(qr)                  
