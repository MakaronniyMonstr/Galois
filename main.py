import sympy as sp
from sympy.abc import x, y
from tabulate import tabulate

import galoislab


def el_mod(p, modulus):
    for idx, _ in enumerate(p):
        p[idx] = int(p[idx] % modulus)
    return p


def create_table(elements, func):
    table = [elements]

    for i in elements:
        row = [i]

        for j in elements:
            row.append(func(i, j))
        table.append(row)

    print(tabulate(table, headers='firstrow', tablefmt='fancy_grid'))


def main():
    irr = x**3+x+1
    p = 2
    m = 3
    elements = [0]

    for i in range(1, p ** m):
        res = sp.rem(x**i, irr, modulus=p, symmetric=False)
        print(res)

        if res not in elements:
            elements.append(res)
        else:
            print('ALARM')

    print(len(elements))

    #create_table(elements.copy(), lambda p1, p2: sp.rem(p1 * p2, irr, modulus=p, symmetric=False))
    #create_table(elements.copy(), lambda p1, p2: sp.rem(p1 + p2, x**p, modulus=p, symmetric=False))


if __name__ == '__main__':

    gf = galoislab.GaloisLab(x=x, p=2, n=3, irr=sp.Poly(x**3+x+1, x))

    gf.extend(y=y, n=2, ext_irr=sp.Poly(y ** 2 + (x + 1) * y + x ** 2 + x + 1, y, x))
    gf.print_ext()
