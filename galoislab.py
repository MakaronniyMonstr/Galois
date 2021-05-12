import sympy as sp
from tabulate import tabulate


class GaloisLab:
    def __init__(self, x: sp.Symbol, p: int, n: int, irr: sp.Poly):
        # Extension parameters
        self.ext_elements = None
        self.ext_irr = None
        self.y = None

        self.p = p
        self.n = n
        self.irr = irr
        self.x = x

        self.dim = p ** n
        # Initial elements
        self.elements = [sp.Poly(0, self.x, modulus=self.p)] + \
                        [sp.Poly(self.x ** i, self.x, modulus=self.p) for i in range(n)]

        for i in range(n, self.dim):
            el = sp.rem(
                sp.Poly(self.x ** i, x), self.irr,
                modulus=p, symmetric=False
            )

            if el not in self.elements:
                self.elements.append(el)
            else:
                if el == 1 and len(self.elements) == self.dim:
                    print('Last element is 1. Irreducible polynomial is really primitive.')
                else:
                    raise ValueError(
                        'Repeated field element detected; please make sure your irreducible polynomial is primitive.')

    def extend(self, y: sp.Symbol, n: int, ext_irr: sp.Poly):
        self.y = y
        self.ext_irr = ext_irr
        self.ext_elements = [sp.Poly(0, self.ext_irr.gens[::-1], modulus=self.p)] + \
                            [sp.Poly(self.y ** i, self.ext_irr.gens[::-1], modulus=self.p) for i in range(n)]

        for i in range(n, self.dim ** n):

            el = sp.rem(
                sp.Poly(self.y ** i, self.y), self.ext_irr,
                modulus=self.p, symmetric=False
            )

            el = sp.rem(
                sp.Poly(el, el.gens[::-1]), self.irr,
                modulus=self.p, symmetric=False
            )

            if el not in self.ext_elements:
                self.ext_elements.append(el)
            else:
                if el == 1 and len(self.elements) == self.dim:
                    print('Last element is 1. Irreducible polynomial is really primitive.')
                else:
                    raise ValueError(
                        'Repeated field element detected; please make sure your irreducible polynomial is primitive.')

    def print(self):
        print(f'{len(self.elements)} elements:')

        self._create_table(self.elements, lambda p1, p2: p1 + p2)
        self._create_table(self.elements, self._mult)

    def print_ext(self):
        print(f'{len(self.ext_elements)} elements:')

        for el in self.ext_elements:
            print(f'{el.args[0]}')

        self._create_table(self.ext_elements, lambda p1, p2: p1 + p2)
        self._create_table(self.ext_elements, self._mult_ext)

    def _mult_ext(self, p1, p2):
        res = p1 * p2

        res = sp.rem(
            sp.Poly(res, res.gens[::-1]), sp.Poly(self.ext_irr),
            modulus=self.p, symmetric=False
        )

        res = sp.rem(
            sp.Poly(res, res.gens[::-1]), sp.Poly(self.irr),
            modulus=self.p, symmetric=False
        )

        return res

    def _mult(self, p1, p2):
        return sp.rem(p1 * p2, self.irr, modulus=self.p, symmetric=False)

    def _create_table(self, elements, func):
        table = [[el.args[0] for el in elements]]

        for i in elements:
            row = [i.args[0]]

            for j in elements:
                row.append(func(i, j).args[0])
            table.append(row)

        print(tabulate(table, headers='firstrow', tablefmt='fancy_grid'))
