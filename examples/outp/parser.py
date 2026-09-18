#%%
import re
import copy
import math
from sympy import symbols, zoo, Rational, Integer
import matplotlib.pyplot as plt
import matplotlib as mpl
from pathlib import Path
from cycler import cycler
mpl.use("Agg")

HARTREE_TO_EV = 27.211386245988


class Parser_OUTP:
    """
    Docstring
    ---------
    a = Parser_OUTP(line=((0,0,0),(1,1,0)), path='mno2afm.outp')
    a.extract_outp()
    a.extract_energy()
    a.save_graph([40, 41], format='pdf')
    """
    def __init__(self, line: tuple, path: str):
        self.line = line
        self.path = path
        self.points_in_line = []

    def extract_points(self, coordinates: list):
        results = []
        for points in coordinates:

            if points:

                result = re.findall(r'\(([\d\s]*)\)', points)
                for j in result:
                    results.append(tuple([int(i) for i in j.strip().split()]))

        return results

    def decide_point(self, points: tuple):
        """
        With functions decide which points belong to definded line.

        Parameters
        ----------
        line : tuple
            Contain coordinats 2 points that define a line.
        points : tuple
            Coordinats of points in kvec.
            """

        x, y, z = symbols("x, y, z")
        xyz = [x, y, z]

        line = self.line

        exprx = (x - line[0][0]) / (line[1][0] - line[0][0])
        expry = (y - line[0][1]) / (line[1][1] - line[0][1])
        exprz = (z - line[0][2]) / (line[1][2] - line[0][2])

        ex = [exprx, expry, exprz]
        ex_bool = [
            True if not i.has(zoo) else False for i in [exprx, expry, exprz]
        ]

        true_points = []

        # TODO Работает только если число точек больше 1, исправить добавив блок try
        for point in points:
            if point == (0, 0, 0) and point in line:
                true_points.append(point)
            elif all(
                [bool(i) == bool(j) for i, j in list(zip(point, ex_bool))]):
                if len(
                        set([
                            ex[exp_number].subs([(xyz[exp_number],
                                                  point[exp_number])])
                            for exp_number in range(len(ex))
                            if ex_bool[exp_number]
                        ])) == 1:
                    true_points.append(point)
            else:
                pass

        del line

        return true_points

    def extract_outp(self):
        """Извлекает данные из файла outp
        """
        n_electron_regex = r'N\. OF ELECTRONS PER CELL\s*([\d]*)\s*[\s+\w+\(+\)+\*-]+'
        # TODO Регулярное выражение работает не верно, не находит фактор, если он меньше 10
        coordinate_factor_regex = r'K POINTS COORDINATES \(OBLIQUE COORDINATES IN UNITS OF IS = ([\d]*)\)'
        # n_points_regex = r'POINTS IN THE IBZ\s*([\d]+)\s'
        _point_regex = r' EIGENVALUES\s[-\s*\w\=]+\(([\d\s]*)\)'

        self.n_electrons, self.n_points, self.coordinate_factor = int(), int(
        ), int()
        k_points_coords = []
        k_point_counter = False

        with open(self.path, 'r') as f:
            for line in f:
                # TODO Добавить if result: ... везде
                if 'N. OF ELECTRONS PER CELL' in line:
                    result = re.search(n_electron_regex, line)
                    self.n_electrons = int(result.group(1))
                elif 'NUMBER OF K POINTS IN THE IBZ' in line:
                # TODO Изменить на поиск с регулярным выражением
                    self.n_points = int(line.split()[7])
                elif 'K POINTS COORDINATES' in line:
                    result = re.search(coordinate_factor_regex, line)
                    self.coordinate_factor = int(result.group(1))
                    # append coord line to list

                    _line = line
                    while _line != '\n':
                        _line = f.readline()
                        k_points_coords.append(_line)

                    k_point_counter = True

                elif k_point_counter:
                    _k_points_coords = tuple(
                        self.extract_points(k_points_coords))
                    k_points_coords = copy.deepcopy(_k_points_coords)
                    del _k_points_coords
                    k_point_counter = False
                    self.points_in_line = self.decide_point(k_points_coords)
                    break

    def points_spliter(self, new_points: str):
        _new_points = []
        lst = [int(i) for i in new_points.strip().split()]
        sz = 3
        if len(lst) % 3 == 0:
            for j in range(0, len(lst), sz):
                _new_points.append(tuple(lst[j:j + sz]))
        else:
            _new_points = [-1]  # TODO Запомнить, проверка на ошибку
        return _new_points

    def add_points(self, new_points: str):
        for i in self.points_spliter(new_points):
            if i in self.points_in_line or i == -1:
                pass
            else:
                self.points_in_line += [i]

    def delete_points(self, del_points):
        for i in self.points_spliter(del_points):
            if i in self.points_in_line:
                self.points_in_line.remove(i)
            else:
                pass

    def split_energy(self, line: str, EVBM: float = 0):
        """Извлекает значение энергии электронных 
        уровней из строки Crystal, переводит значения в eV 
        и вычитает из значений EVBM.

        Args:
            line: str 
                Строка содержащая значения энергии для электронных уровней
            EVBM: float
                Максимальное значение верхнего электронного уровня в eV

        Returns:
            list: Список с энергниями
        """
        _regex = r'(-?\d*\.\d*E[-+]+\d{2})'
        results = []

        if line:
            result = re.findall(_regex, line)
            for j in result:
                results.append((float(j.strip()) * HARTREE_TO_EV) - EVBM)

        return results

    def extract_energy(self, EVBM: float = 0):
        """Вытаскивает энергию отдельных точек
        """
        _point_regex = r' EIGENVALUES\s[-\s*\w\=]+\(([\d\s]*)\)'
        self.alpha, self.beta = {}, {}
        alpha_switch = False
        beta_switch = False
        point_switch = False

        with open(self.path, 'r') as f:
            for line in f:
                if '    ALPHA      ELECTRONS' in line:
                    alpha_switch = True
                elif '    BETA       ELECTRONS' in line:
                    alpha_switch = False
                    beta_switch = True

                elif 'EIGENVALUES' in line:
                    _point = re.search(_point_regex, line)
                    _point = tuple(
                        [int(i) for i in _point.group(1).strip().split()])
                    if _point in self.points_in_line:
                        point_switch = True
                        _ = []
                        if alpha_switch:
                            self.alpha[_point] = _
                        elif beta_switch:
                            self.beta[_point] = _
                        else:
                            continue
                    else:
                        point_switch = False

                elif point_switch and line != '\n' and (alpha_switch
                                                        or beta_switch):
                    _ += self.split_energy(line, EVBM)

    def form_datasets(
        self,
        e_indexes: tuple,
    ):
        """Подготавливает данные для отрисовки 
        в matplotlib формируя массивы данных по х и у

        Args:
            e_indexes (tuple):
                Набор уровней, значения которых нужно вытащить из словоря

        Returns:
            Список с расстоянием от начальнгой точки до 
            текущей и два списка с значениями энергий в электронвольтах
        """
        self.points_in_line.sort()
        data_x = [
            self.get_length(point_1=self.points_in_line[0],
                            point_2=self.points_in_line[i])
            for i in range(len(self.points_in_line))
        ]
        data_y_a = []
        data_y_b = []

        for index in e_indexes:
            _a = [self.alpha[i][index - 1] for i in self.points_in_line]
            data_y_a.append(_a)
            _b = [self.beta[i][index - 1] for i in self.points_in_line]
            data_y_b.append(_b)

        return data_x, data_y_a, data_y_b

    def get_length(self, point_1: tuple, point_2: tuple):
        return (math.sqrt(
            sum([(point_1[i] - point_2[i])**2 for i in range(len(point_1))])))

    def sym_points(self):
        sym_points = []
        for i in self.points_in_line:
            _ = []
            for j in i:
                if j == 0:
                    _.append(Integer(j))
                else:
                    _.append(Rational(j, self.coordinate_factor))
            sym_points.append(tuple(_))
        return sym_points

    def save_graph(self,
                   e_indexes: tuple,
                   format: str = None,
                   interpolation: bool = False,
                   figsize: str = None):

        import numpy as np
        from scipy import interpolate
        from scipy import optimize

        if format:
            format = format
        else:
            format = 'png'

        if figsize and len(figsize.strip().split()) == 2:
            try:
                figsize = [int(i) for i in (figsize.strip().split())]
            except:
                figsize = (7, 6)
        else:
            figsize = (7, 6)

        default_cycler = (
            cycler(color=['red', 'blue'] * 6) + cycler(linestyle=[
                '-', '-', '--', '--', (0, (3, 1, 1, 1)),
                (0, (3, 1, 1, 1)), ':', ':', '-.', '-.', (0,
                                                          (1, 1)), (0, (1, 1))
            ]))
        plt.rc('axes', prop_cycle=default_cycler)

        data_x, data_y_a, data_y_b = self.form_datasets(e_indexes)
        fig = plt.figure(figsize=figsize)
        ax = fig.add_axes([0.1, 0.1, 0.8, 0.8])  # main axes

        # Обратный порядорк чтобы в легенде первыми шли уровни с высоким номером
        for i in range(len(e_indexes) - 1, -1, -1):
            if interpolation == True:
                f_alpha = interpolate.interp1d(data_x, data_y_a[i], 'cubic')
                f_beta = interpolate.interp1d(data_x, data_y_b[i], 'cubic')
                new_x = np.linspace(min(data_x), max(data_x), 200)
                #Line
                ax.plot(new_x, f_alpha(new_x), label=f'{e_indexes[i]} (Up)', linewidth=0.4)
                ax.plot(new_x, f_beta(new_x), label=f'{e_indexes[i]} (Down)', linewidth=0.4)
                #Dots
                ax.scatter(data_x, data_y_a[i], marker='o', s=0.5)
                ax.scatter(data_x, data_y_b[i], marker='o', s=0.5)
            else:
                ax.plot(data_x, data_y_a[i], label=f'{e_indexes[i]} (Up)')
                ax.plot(data_x, data_y_b[i], label=f'{e_indexes[i]} (Down)')

        ax.set_xticklabels(self.sym_points()[0::len(self.sym_points()) - 1],
                           fontsize=12)
        ax.tick_params(axis='y', labelsize=12)
        ax.tick_params(direction="in")
        ax.set_xlim(xmin=min(data_x), xmax=max(data_x))
        ax.set_xticks(data_x[0::len(self.sym_points()) - 1])
        ax.set_xlabel('Wave Vector', fontsize=16)
        ax.set_ylabel('E-E$_{VBM}$, eV', fontsize=16)
        ax.legend()

        plt.savefig(fname=Path(self.path).with_suffix(f'.{format}'),
                    format=format,
                    bbox_inches='tight',
                    pad_inches=0,
                    dpi=360,
                    transperent=True)
        
        plt.close(fig)

    def save_txt(self, e_indexes: tuple):
        data_x, data_y_a, data_y_b = self.form_datasets(e_indexes)

        with open(Path(self.path).with_suffix(f'.txt'), 'w') as f:
            f.write(f'Point/Index Alpha Beta\n')
            for j in range(len(e_indexes)):
                f.write(f'{e_indexes[j]}\n')
                for i in range(len(data_y_a[j])):
                    f.write(
                        f'{self.points_in_line[i]} {data_y_a[j][i]} {data_y_b[j][i]}\n'
                    )


# %%
