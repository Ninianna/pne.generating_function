from cgi import test
from contextlib import nullcontext
from os import read
import polytope.polytope as pol
import numpy as np
import fractions as fra
import os
import regex as re
from Schubert3 import *
import sympy as sp
import random
import math
import pandas as pd
from fractions import Fraction
import time

class GF:
    def __init__(self,sign,numers,denoms):
        self.sign = sign
        self.numers = numers
        self.denoms = denoms

def gbs(s):
    a,b = s[0],s[1]
    a = a // math.gcd(a, b) * b // math.gcd(a, b) * math.gcd(a, b)
    if len(s)>2:
        for i in range(2,len(s)):
            b = s[i]
            a = a//math.gcd(a,b) * b//math.gcd(a,b) * math.gcd(a, b)
    return a

def read_to_latte(file_name,gf_file_name):
    print("read file")
    with open(file_name) as f:
        lines = f.read().splitlines()
    A, b = [], []
    file = open(gf_file_name,'w')
    for l in lines:
            l = l.split()
            A.append([int(s) for s in l[:-1]])
            b.append(int(l[-1]))
    file.write(str(len(A)) + " " + str(len(A[0])+1) + "\n")
    for i,v in enumerate(b):
        file.write(str(v) + " ")
        for j,m in enumerate(A[i]):
            A[i][j] = m*(-1)
            file.write(str(A[i][j]) + " ")
        file.write("\n")
    file.close()
    rename_output = gf_file_name.replace('.txt','.hrep.latte')
    new_file = os.rename(gf_file_name,rename_output)
    
def read_polytope(file_name):
    with open(file_name) as f:
        lines = f.read().splitlines()
    
    A, b = [], []
    for l in lines:
            l = l.split()
            A.append([int(s) for s in l[:-1]])
            b.append(int(l[-1]))
    polytope = pol.Polytope(np.array(A),np.array(b),normalize=False)
    print("finished reading polytope")
    return polytope

def project_polytope(poly,dim):
    proj_poly = pol.projection(poly,dim)
    print("finished projection")
    ieq = []
    for i,v in enumerate(proj_poly.b):
        A = []
        A.append(Fraction(str(round(v,2))).limit_denominator())
        for j,m in enumerate(proj_poly.A[i]):
            temp = m*(-1)
            A.append(Fraction(str(round(temp,2))).limit_denominator())
        ieq.append(A)   
    return ieq

def write_to_latte(ieq, file_name):
    print("starting write to Latte")
    file = open(file_name,'w')
    poly = Polyhedron(vertices=ieq, base_ring=QQ)
    poly.base_ring()
    h_poly = poly.cdd_Hrepresentation()
    h = h_poly.splitlines()
    file.write(str(h[2][1:-9])+"\n")
    h = h[3:-1]
    for i in h:
        for j,v in enumerate(i):
            if(j==0):
                continue
            else:
                file.write(v)
        file.write("\n")
    file.close()
    rename_output = file_name.replace('.txt','.hrep.latte')
    new_file = os.rename(file_name,rename_output)
    print("finish write to Latte")

# please change your own path
def generate_function(file_name):
    _str_1 = "cp -r YourPath"+ file_name + " ~/YourLattePath;"
    os.system(_str_1)
    _str_3 = "cd YourLattePath; ./count --multivariate-generating-function " + file_name + ";"
    os.system(_str_3)
    _str_4 = "cp -r YourLattePath" + file_name + ".rat /YourPath;"
    os.system(_str_4)
    new_file = file_name + ".rat"
    rename_output = new_file.replace('.hrep.latte.rat','.txt')
    os.rename(new_file, rename_output)

# please change your own path
def generate_unbounded(file_name):
    _str_1 = "cp -r YourPath"+ file_name + " ~/YourLattePath;"
    os.system(_str_1)
    _str_3 = "cd /YourLattePath; count --compute-vertex-cones=4ti2 --multivariate-generating-function " + file_name + ";"
    os.system(_str_3)
    _str_4 = "cp -r /YourLattePath" + file_name + ".rat /YourPath;"
    os.system(_str_4)
    new_file = file_name + ".rat"
    rename_output = new_file.replace('.hrep.latte.rat','.txt')
    os.rename(new_file, rename_output)

def Enu(file_name, generate_file_name):
    # input file
    time_start = time.time()
    with open(file_name) as f:
        lines = f.read().splitlines()
    
    A, b = [], []
    for l in lines:
            l = l.split()
            A.append([int(s) for s in l[:-1]])
            b.append(int(l[-1]))
    for i,v in enumerate(b):
        for j,m in enumerate(A[i]):
            A[i][j] = m*(-1)
        A[i].insert(0,b[i])
    poly = Polyhedron(ieqs=A, base_ring=QQ)
    print("starting enumeration")
    enu = poly.integral_points()
    print("finish enumeration")
    file = open(generate_file_name, 'w')
    enu_new = []
    for i in enu:
        i = list(i)
        enu_new.append(i)
        file.write(str(i)+'\n')
    file.close()
    time_end = time.time()
    t = time_end - time_start
    print("Enu:",t)
    return enu_new, generate_file_name

def Proj(point, proj_vector, generate_file_name):
    print("starting projection")

    file = open(generate_file_name, 'w')
    points = []
    for i in point:
        temp = list(i)
        point = []
        for j,v in enumerate(proj_vector):
            if(v == 1):
                point.append(temp[j])
        file.write(str(point)+'\n')
        points.append(point)
    file.close()
    final_point = []
    for i in points:
        if(i not in final_point):
            final_point.append(i)
    return final_point,generate_file_name

def Union(vectors, generate_file_name):
    F = []
    for i in vectors:
        for j in i:
            if(j not in F):
                F.append(j)

    print("The points of union: ",len(F))
    file = open(generate_file_name, 'w')
    for i in F:
        file.write(str(i)+'\n')
    file.close()
    return F,generate_file_name

# please change your own path
def generate_unbounded_function(file_name,N):
    _str_1 = "cp -r YourPath"+ file_name + " ~/YourLattePath;"
    os.system(_str_1)
    _str_3 = "cd /YourLattePath; count --compute-vertex-cones=4ti2 --multivariate-generating-function " + file_name + ";"
    os.system(_str_3)
    _str_4 = "cp -r /YourLattePath" + file_name + ".rat YourPath;"
    os.system(_str_4)
    new_file = file_name + ".rat"
    rename_output = new_file.replace('.hrep.latte.rat','.txt')
    os.rename(new_file, rename_output)
    generate_fun = file_format("test.txt")
    class_gf = generate_class(generate_fun,N)

    return class_gf

def file_format(file_name):
    with open(file_name) as f:
        lines = f.read().splitlines()
    f_ = []
    positive, negative = [], []
    for i,v in enumerate(lines):
        if(i==0):
            if(v[:5] == '(-1)*'):
                negative.append(v[5:])
            else:
                positive.append(v)
            continue
        else:
            if(v[3:8] == '(-1)*'):
                negative.append(v[8:])
            else:
                positive.append(v[3:])
    return positive, negative

def split_generating_function(gen_fun):
    formula_pair = []
    for formula in gen_fun:
        depth = 0
        for i in range(len(formula)-1):
            if(formula[i] == "("):
                depth += 1
            elif(formula[i]==")"):
                depth -=1
            elif(depth == 0 and formula[i]=="/"):
                formula_pair.append((formula[:i] ,formula[i+1:]))
    return formula_pair

def get_vectors(gen_fun_pair,N):
    vectors = []
    for i in gen_fun_pair:
        denom = i[1]
        v = re.findall(r'(1-(x\[(\d+)\](\^\(([\+\-]?\d+)\))?\*?)+)', denom)
        vector = []
        for j in v:
            vec = [0 for _ in range(N)]
            vec_ = re.findall(r'(x\[(\d+)\](\^\(([\+\-]?\d+)\))?)', j[0])
            for d in vec_:
                vec[int(d[1])] = 1 if d[3]=='' else int(d[3])
            vector.append(vec)
        vectors.append(vector)
    return vectors

def get_vertices(gen_fun_pair,N):
    vertices = []
    for i in gen_fun_pair:
        numer = i[0]
        vertice = [0 for _ in range(N)]
        ver = re.findall(r'(x\[(\d+)\](\^\(([\+\-]?\d+)\))?)', numer)
        for d in ver:
            vertice[int(d[1])] = 1 if d[3]=='' else int(d[3])
        vertices.append(vertice)
    return vertices

def generate_class(generate_fun,N):
    class_gf = []
    positive = split_generating_function(generate_fun[0])
    negative = split_generating_function(generate_fun[1])
    if(len(positive)!=0):
        for i in range(len(positive)):
            temp = GF(1,get_vertices(positive,N)[i],get_vectors(positive,N)[i])
            class_gf.append(temp)
    if(len(negative)!=0):
        for j in range(len(negative)):
            temp = GF(-1,get_vertices(negative,N)[j],get_vectors(negative,N)[j])
            class_gf.append(temp)
    return class_gf

def enumerate_points(vertices, file_name):
    poly = Polyhedron(vertices=vertices, base_ring=QQ)
    poly.base_ring()
    enu = poly.integral_points()
    points = []
    for i in enu:
        if(list(i) not in vertices):
            points.append(i)
    file = open(file_name, 'w')
    for j in points:
        file.write(str(j)+"\n")
    file.close()
    return points

def negative_to_positive(vertices, vectors):
    for i,v in enumerate(vectors):
        v[0] = [x*(-1) for x in v[0]]
        vertices[i] = [vertices[i][j] + v[0][j] for j in range(len(vertices[i]))]
    return vertices

def union(gf_p, gf_inter):
    positive = None
    negative = None
    for i in gf_p:
        g = file_format(i)
        positive += g[0]
        negative += g[1]
    positive += gf_inter[1]
    negative += gf_inter[0]
    new_positive = []
    new_negative = []
    for i in positive:
        if i not in negative:
            new_positive.append(i)
    for i in negative:
        if i not in positive:
            new_negative.append(i)

    return new_positive, new_negative

import pickle

def substract(file_1, file_2, d):
    f1 = file_format(file_1)
    f2 = file_format(file_2)
    positive = f1[0] + f2[1]
    negative = f1[1] + f2[0]
    new_positive, new_negative = [], []
    for i in positive:
        if i not in negative:
            new_positive.append(i)
    for i in negative:
        if i not in positive:
            new_negative.append(i)
    # print(negative,positive)
    gen_fun_1 = split_generating_function(new_positive)
    gen_fun_2 = split_generating_function(new_negative)
    generating_function = []
    for i in range(len(new_positive)):
        temp = GF(1,get_vertices(gen_fun_1,d)[i],get_vectors(gen_fun_1,d)[i])
        generating_function.append(temp)
    for j in range(len(new_negative)):
        temp = GF(-1,get_vertices(gen_fun_2,d)[j],get_vectors(gen_fun_2,d)[j])
        generating_function.append(temp)
    return generating_function

def intersection(ieqs):
    for i in ieqs:
        for j in i:
            if(j not in ieqs[0]):
                ieqs[0].append(j)
    return ieqs[0]

def get_union_vertice(pos_gf, neg_gf, d):
    pos_ver = get_vertices(split_generating_function(pos_gf), d)
    neg_ver = get_vertices(split_generating_function(neg_gf), d)
    neg_vec = get_vectors(split_generating_function(neg_gf), d)
    neg_ver2 = negative_to_positive(neg_ver, neg_vec)
    new_ver = pos_ver + neg_ver2
    return new_ver

def Todd_Class(variable, variable_list, dimension):
    R = PolynomialRing(QQ,variable,order=TermOrder('wdegrevlex',variable_list))
    f = 1
    for i in R.gens():
        f = f + i
    print(f)
    g = logg(f,dimension)+dimension
    print("g=",g)
    t_class = todd(g,dimension)
    return t_class

def find_vector(h,g,dimension):
    k = dimension
    l = []
    b1 = []
    b2 = []
    for _h in h:
        for denom in _h.denoms:
            b1.append(denom)
    for _g in g:
        for denom in _g.denoms:
            b2.append(denom)

    while(True):
        if(len(l) != 0):
            break
        else:
        #create random l
            dot_value = []
            for _ in range(k):
                l.append(random.randint(-9,9))
            #check l with b
            for _b1 in b1:
                bi = np.array(_b1).reshape(len(_b1),1)
                biT = bi.transpose()
                value1 = float(np.dot(biT,l))
                dot_value.append(value1)
            for _b2 in b2:
                bi = np.array(_b2).reshape(len(_b2),1)
                biT = bi.transpose()
                value2 = float(np.dot(biT,l))
                dot_value.append(value2)
            for i in dot_value:
                if(i == 0 or i>=4 or i<=-4):
                    l = []
                    break
    print(l)

def Residue(h,variable_string,variable_list,dimension):
    k = dimension

    l = []
    b = []
    for _h in h:
        for denom in _h.denoms:
            if(denom != [0,0,0,0]):
                b.append(denom)
    value_ = []
    while(True):
        if(len(l) != 0):
            break
        else:
        #create random l
            for _ in range(k):
                l.append(random.randint(-11,11))
            #check l with b
            for _b in b:
                bi = np.array(_b).reshape(len(_b),1)
                biT = bi.transpose()
                v_ = float(np.dot(biT,l))
                value_.append(abs(v_))
                if(v_ == 0):
                    # or value>=4 or value<=-4
                    value_ = []
                    l = []
                    break    
    #Find e
    e = []
    for _h in h:
        ei = []
        if(len(_h.denoms)==dimension):
            for denom in _h.denoms:
                bi = np.array(denom).reshape(len(denom),1)
                biT = bi.transpose()
                ei.append(int(np.dot(biT,l)))
            e.append(ei)
        else:
            for denom in _h.denoms:
                bi = np.array(denom).reshape(len(denom),1)
                biT = bi.transpose()
                ei.append(int(np.dot(biT,l)))
                for j in range(dimension - len(_h.denoms)):
                    temp = 0
                    ei.append(temp)
            e.append(ei)
    #Find t
    t = []
    for h_ in h:
        ti = []
        # if(len(_h.denoms)==dimension):
        numer = h_.numers
        ai = np.array(numer).reshape(len(numer),1)
        aiT = ai.transpose()
        value = int(np.dot(aiT,l))
        value = value*(_h.sign)
        ti.append(value)
        t.append(ti)

    print(e,t)
    print(len(e),len(t))
    variable = variable_string
    td_c = []
    for i in range(k):
        if(i==0):
            td =  Todd_Class(variable[:k*2],variable_list[:k],k) - Todd_Class(variable[:(k*2)-(1+(2*(i+1)))],variable_list[:k-(i+1)],k-(i+1))
            td_c.append(td)
        else:
            if(k-(i+1)!=0):
                td =  Todd_Class(variable[:(k*2)-(1+(2*i))],variable_list[:k-i],k-i) - Todd_Class(variable[:(k*2)-(1+(2*(i+1)))],variable_list[:k-(i+1)],k-(i+1))
                td_c.append(td)
            else:
                td =  Todd_Class(variable[:(k*2)-(1+(2*i))],variable_list[:k-i],k-i)-1
                td_c.append(td)
    td_c.append(1)
    E = []
    TD = []
    Ans = 0
    for i,p in enumerate(e):
        num = 1
        for q in p:
            num = num*q
        if(num != 0):
            count = h[i].sign / num
            E.append(count)
        else:
            count = 0
            E.append(count)

        temp = 0
        for n in range(k+1):
            if(n!=k):
                z = td_c[n].variables()
                d = len(z)
                if(len(p)>=d):
                    z = tuple(p[:k-n])
                else:
                    redunt = []
                    for s in p:
                        redunt.append(s)
                    for i in range(d-len(p)):
                        redunt.append(0)
                    z = tuple(redunt)
                temp = temp + (pow(t[i][0],n)/factorial(n))*td_c[n](z)
            else:
                temp = temp + (pow(t[i][0],n)/factorial(n))
        TD.append(temp)
        Ans = Ans + count*temp
    return str(l),float(Ans)

def generate(upper_bound, lower_bound, dimension, var_index):
    file_name = "test.hrep.latte"
    file = open(file_name,'w')
    file.write(str(2+dimension-1)+" "+str(dimension+1)+"\n")
    for d in range(2):
        for i in range(dimension+1):
            if(i == 0):
                if(d == 0):
                    file.write(str(upper_bound)+" ")
                else:
                    file.write(str(lower_bound)+" ")
            elif(i == var_index):
                if(d == 0):
                    file.write("-1 ")
                else:
                    file.write("1 ")
            else:
                file.write("0 ")
        file.write("\n")
    for d in range(1, dimension+1):
        line = [0 for _ in range(dimension+1)]
        if(d != var_index):
            line[d] = 1
            for j in line:
                file.write(str(j)+" ")
            file.write("\n")
    file.close()
    return file_name
  
def condition(v1, v2, l):
    for i,v in enumerate(v1.denoms):
        if(np.dot(l, v) > 0):
            v1.denoms[i] = [i*(-1) for i in v]
            v1.numers = [v1.numers[j] + v1.denoms[i][j] for j in range(len(v1.denoms[i]))]
            v1.sign = v1.sign * -1  # 乘上負號
    for i,v in enumerate(v2.denoms):
        if(np.dot(l, v) > 0):
            v2.denoms[i] = [i*(-1) for i in v]
            v2.numers = [v2.numers[j] + v2.denoms[i][j] for j in range(len(v2.denoms[i]))]
            v2.sign = v2.sign * -1  # 乘上負號
    return v1, v2

def product(m1, m2):
    # print(m1.numers, m1.denoms)
    # print(m2.numers, m2.denoms)
    Y = []  # <=0
    Z = []  # =0
    H = []
    T = []

    # 非負限制式
    for n in range(0, len(m1.denoms)+len(m2.denoms)):
        z = []
        z.append(0)
        for m in range(0, len(m1.denoms)+len(m2.denoms)):
            if n == m:
                z.append(1)
            else:
                z.append(0)
        z_str = " ".join(map(str, z))
        Y.append(z_str)

    # =0,[b-A]
    for i in range(len(m1.numers)):
        x = []
        b = m2.numers[i]-m1.numers[i]
        x.append(b)
        for j in range(len(m1.denoms)):
            temp1 = m1.denoms[j][i] * -1
            x.append(temp1)
        for k in range(len(m2.denoms)):
            temp2 = m2.denoms[k][i]
            x.append(temp2)
        x_str = " ".join(map(str, x))
        Z.append(x_str)
        #temp = x_arr * -1
        # Y.append(temp)

    h = []
    h.append(len(Y)+len(Z))
    h.append(len(m1.denoms)+len(m2.denoms)+1)
    h_str = " ".join(map(str, h))
    H.append(h_str)

    t = ["linearity"]
    t.append(len(Z))
    for i in range(0, len(Z)):
        t.append(len(Y) + 1 + i)
    t_str = " ".join(map(str, t))
    T.append(t_str)

    return H, Y, Z, T

def substitution(v_1, arr, sign):
    print(v_1.numers, v_1.denoms)
    dim = len(v_1.numers)
    if(len(v_1.denoms) < dim):
        for d in range(dim-len(v_1.denoms)):
            vec = [0 for _ in range(dim)]
            v_1.denoms.append(vec)
    class_gf = []
    for i in arr:
        #numers
        N = [ 0 for i in range(dim)]
        for j,v in enumerate(i.numers[:dim]):
            for k,n in enumerate(v_1.denoms[j]):
                numer = v*n
                N[k] += numer
        for m in range(dim):
            N[m] +=  v_1.numers[m]
        #denoms
        D = []
        for j in i.denoms:
            denoms = [0 for i in range(dim)]
            for k,v in enumerate(j[:dim]):
                for n,d in enumerate(v_1.denoms[k]):
                    denom = v*d
                    denoms[n] += denom
            D.append(denoms)
        gf = GF(i.sign*sign,N,D)
        class_gf.append(gf)
    return class_gf 

def hadamard(p1, p2, l):
    con_p1, con_p2 = condition(p1, p2, l)
    sign = con_p1.sign * con_p2.sign
    H, Y, Z, T = product(con_p1, con_p2)
    txt = H + Y + Z + T
    fo = open("test.hrep.latte", "w")
    fo.writelines([line+'\n' for line in txt])
    fo.close()
    N = len(con_p1.numers) + len(con_p2.numers)
    
    new_generating = generate_unbounded_function("test.hrep.latte", N)
    print("new:",new_generating)
    if (len(new_generating)==0):
        print("None")
        result = []
    else:
        result = substitution(con_p1, new_generating, sign)
    return result

def Hadamard(a,b):
    F,l,p1,p2 = [],[],[],[]
    for i in a:
        for _i in i.denoms:
            p1.append(_i)
    for j in b:
        for _j in j.denoms:
            p2.append(_j)
    while(True):
        if(len(l) != 0):
            break
        else:
            for i in range(len(a[0].numers)):
                l.append(random.randint(1, 10))
            for j in p1:
                value1 = np.dot(l,j)
            for k in p2:
                value2 = np.dot(l,k)
            if(value1 == 0 and value2 == 0):
                l = []
    print(l)
    # for i in range(len(a)):
    for i in range(1):
        # for j in range(len(b)):
        for j in range(1):
            h = hadamard(a[i], b[j], l)
            if(h!=[]):
                for k in h:
                    F.append(k)
    print("FFFF:",F)
    return F

def Finding_actual_bound(upper_bound, reference_bound, rational_generating_fuction,dimension,var_index,var,var_list):
    u = upper_bound
    r = reference_bound
    f = rational_generating_fuction
    while True:
        _reference = int(math.ceil((u + r) / 2))
        if(_reference == u):
            break
        h = Hadamard(generate_unbounded_function(generate(_reference, u,dimension,var_index),dimension), f)
        for i in h:
            print(i.sign,i.numers,i.denoms)
        R = Residue(h,var,var_list,dimension)
        if( R[1] < 0):
            u = _reference
        else:
            if(R[1] >= 0 ):
                r = _reference
            else:
                u = r
    return u

def cut_bounds(upper, refer, lower, rgf, dim, var, var_list):
    ubounds = []
    lbounds = []
    for i in dim:
        ubound = Finding_actual_bound(upper,refer,rgf,dim,i+1,var,var_list)
        lbound = Finding_actual_bound(refer,lower,rgf,dim,i+1,var,var_list)
        ubounds.append(ubound)
        lbounds.append(lbound)
    return ubounds, lbounds

def Enumeration(cut_bounds):
    s = cut_bounds[1]
    e = cut_bounds[0]
    def get_point(dim):
        result = []
        if(dim == len(s)):
            result = [""]
        else:
            if(e[dim]<s[dim]):
                for i in range(e[dim],s[dim]+1):
                    for r in get_point(dim+1):
                        result.append(str(i)+","+r)
            else:
                for i in range(s[dim],e[dim]+1):
                    for r in get_point(dim+1):
                        result.append(str(i)+","+r)
        return result
    g = get_point(0)
    exp = []
    for i in g:
        temp = i.split(',')
        v = []
        for j in temp[:-1]:
            v.append(int(j))
        exp.append(v)
    return exp

def NashEquil(vertice, filename, generating_function, var, var_list):
    dimension = len(vertice)
    file_name = filename.replace(".txt",".hrep.latte")
    file = open(file_name,'w')
    file.write(str(dimension*2)+" "+str(dimension+1)+"\n")
    for j in range(2):
        for i,v in enumerate(vertice):
            file.write(str(v) + " ")
            for k in range(dimension):
                if(k == i):
                    if(j==0):
                        file.write("1 ")
                    else:
                        file.write("-1 ")
                else:
                    file.write("0 ")
            file.write("\n")
    file.close()
    generate_unbounded(file_name)
    gf = generate_class(file_format(filename),dimension)
    # for i in gf:
    #     print(i.sign,i.numers,i.denoms)
    h = Hadamard(gf, generating_function)
    file_1 = open("vertice.hrep.latte","w")
    file_1.close()
    file_2 = open("vertice.hrep.latte.rat","w")
    file_2.close()
    file_3 = open("vertice.txt","w")
    file_3.close()
    return h

def validate(players,all,points):
    ne = []
    for i in points:
        if(i not in players):
            if(i in all):
                ne.append(i)
    for i in all:
        if(i not in players):
            if(i not in points):
                ne.append(i)
    return ne

def validation(points,players,all):
    ne = []
    for p in points:
        for i in players:
            cons_status = None
            player_status = None
            player_check = ["" for _ in range(len(players))]
            for j,h in enumerate(i.A):
                for idx, k in enumerate(h):
                    if(k*p[idx] == i.b[j]):
                        break
                else:
                    cons_status = "True"
                if(cons_status != "True"):
                    break
            else:
                player_status = "True"
            if(player_status == "True"):
                player_check[i] = "True"
            else:
                player_check[i] = "False"
        if("False" not in player_check):
            all_status = None
            for m,n in enumerate(all.A):
                for c,d in n:
                    if(d*p[c] != all.b[m]):
                        break
                else:
                    all_status = "True"
                if(all_status != "True"):
                    break
            else:
                ne.append(p)
    return ne              

def txt_to_list(file_name):
    with open(file_name) as f:
        lines = f.read().splitlines()
    vecs = []
    for i in lines:
        temp = i.split('[')
        temp_1 = temp[1].split(']')
        temp_2 = temp_1[0].split(',')
        vec = []
        for j in temp_2:
            vec.append(int(j))
        vecs.append(vec)
    return vecs

def check_duplicate(ver):
    check = []
    for i in ver:
        if(i not in check):
            check.append(i)
    return check


