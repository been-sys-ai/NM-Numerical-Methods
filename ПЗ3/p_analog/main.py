import math
import random
import re
import sys


# Опишем класс реализующий интерфейс сглаживающего сплайна

class Smoothing_spline:

    # параметр сглаживания
    SMOOTH: float
    # точки сетки
    x: list
    # коэффициенты разложения по базису
    alpha: list
    
    a: list
    b: list
    c: list

    # перход на мастер элемнент [-1,1]
    def Trans_To_Mas_El(self, s_num, X):
        Ksi = (
            2.0 * (X - self.x[s_num])
            / (self.x[s_num + 1] - self.x[s_num])
            - 1.0
        )
        return Ksi

    # базисные функции на [-1,1]
    def Bas_Func(self, N, Ksi):
        if N == 1:
            return 0.5 * (1 - Ksi)

        elif N == 2:
            return 0.5 * (1 + Ksi)

        else:
            raise ValueError("Укажите корректный номер функции")

    # производные базисных функций на [-1,1]
    def Der_Bas_Func(self, N, Ksi):
        if N == 1:
            return -0.5

        elif N == 2:
            return 0.5

        else:
            raise ValueError("Укажите корректный номер функции")

    # конструктор сглаживающего сплайна
    def __init__(self, SMOOTH):
        self.SMOOTH = SMOOTH
        self.x = []
        self.alpha = []
        self.a = []
        self.b = []
        self.c = []

    # Функция для обновления сплайна
    def update_spl(self, x_upd, f_d):

        # обновляем точки сетки
        self.x.clear()

        for elem in x_upd:
            self.x.append(elem)

        # определяем число отрезков разбиения
        n = len(self.x) - 1

        # формируем вектор коэффициентов разложения по базису
        self.alpha = self.alpha[:n + 1]
        self.alpha.extend([0.0] * (n + 1 - len(self.alpha)))

        # описываем диагонали матрицы СЛАУ

        self.a = self.a[:n + 1]
        self.a.extend([0.0] * (n + 1 - len(self.a)))

        self.b = self.b[:n + 1]
        self.b.extend([0.0] * (n + 1 - len(self.b)))

        self.c = self.c[:n + 1]
        self.c.extend([0.0] * (n + 1 - len(self.c)))

        # процедура ассемблирования СЛАУ
        def Assembly(i, X, F_Val, w):

            # переход на мастер элемент
            Ksi = self.Trans_To_Mas_El(i, X)
            # вычисление значений баззисных функций на мастер элементе
            f1 = self.Bas_Func(1, Ksi)
            f2 = self.Bas_Func(2, Ksi)
            # внесение вкладов в матрицу СЛАУ
            self.b[i] += (1.0 - self.SMOOTH) * w * f1 * f1
            self.b[i + 1] += (1.0 - self.SMOOTH) * w * f2 * f2

            self.a[i + 1] += (1.0 - self.SMOOTH) * w * f1 * f2
            self.c[i] += (1.0 - self.SMOOTH) * w * f2 * f1
            # внесение вкладов в вектор правой части СЛАУ
            self.alpha[i] += (1.0 - self.SMOOTH) * w * f1 * F_Val
            self.alpha[i + 1] += (1.0 - self.SMOOTH) * w * f2 * F_Val

        # сборка СЛАУ
        for i in range(n):

            W_left = 0.5 if i == 8 else 1.0
            W_right = 0.5 if i + 1 == 8 else 1.0

            Assembly(i, self.x[i], f_d[i], W_left)
            Assembly(i, self.x[i + 1], f_d[i + 1], W_right)

            h = self.x[i + 1] - self.x[i]
            self.b[i] += 1.0 / h * self.SMOOTH
            self.b[i + 1] += 1.0 / h * self.SMOOTH
            self.a[i + 1] -= 1.0 / h * self.SMOOTH
            self.c[i] -= 1.0 / h * self.SMOOTH

        # метод прогонки: прямой ход
        for j in range(1, n + 1):
            self.b[j] -= self.a[j] / self.b[j - 1] * self.c[j - 1]
            self.alpha[j] -= (
                self.a[j] / self.b[j - 1] * self.alpha[j - 1]
            )
        # метод прогонки: обратный ход
        self.alpha[n] /= self.b[n]
        for j in range(n - 1, -1, -1):
            self.alpha[j] = (
                self.alpha[j] - self.alpha[j + 1] * self.c[j]
            ) / self.b[j]

    # Функция для вычисления значения сплайна и его производной в точке
    def get_spl_value_and_dev(self, x_in):

        # машинный ноль
        eps = 1e-7

        # определяем число отрезков разбиения
        n = len(self.x) - 1

        # вектор результата
        res = [0.0] * 3

        for i in range(n):
            if (
                (x_in > self.x[i] and x_in < self.x[i + 1])
                or abs(x_in - self.x[i]) < eps
                or abs(x_in - self.x[i + 1]) < eps
            ):
                # длина отрезка
                local = self.x[i + 1] - self.x[i]
                # переход на мастер элемент
                Ksi = self.Trans_To_Mas_El(i, x_in)
                # вычисление значения сплайна и его производной по базисным функциям
                res[0] = (
                    self.alpha[i] * self.Bas_Func(1, Ksi)
                    + self.alpha[i + 1] * self.Bas_Func(2, Ksi)
                )
                res[1] = (
                    self.alpha[i] * self.Der_Bas_Func(1, Ksi)
                    + self.alpha[i + 1] * self.Der_Bas_Func(2, Ksi)
                ) * 2.0 / local
                res[2] = 0.0

                return res

        raise ValueError("Указано не корректное число")


# Генерация случайных чисел из нормального распределения
def normDen(M, s, N):
    if not math.isfinite(M):
        raise ValueError(
            "Мат ожидание должно быть конечным числом"
        )
    if not math.isfinite(s) or s <= 0.0:
        raise ValueError(
            "Стандартное отклонение должно быть конечным и положительным"
        )

    gen = random
    dist = gen.gauss

    val = []

    for i in range(N):
        val.append(dist(M, s))

    return val


def main(argc, argv):

    # параметры нормального распредления и число точек

    N = 1443
    M = 0.94
    s = 4.28

    # формируем последовательность случайных чисел по нормальному закону

    # блок записи, для возможности проведения анализа весов скалярного произвдеения
    # den = normDen(M, s, N)

    # try:
    #     file = open(
    #         "/Users/bogdanvolkhin/Documents/Учеба ВУЗ/НГТУ/3 курс/ЧМ/ПЗ3/plot/den.txt",
    #         "w",
    #         encoding="utf-8",
    #     )
    # except OSError:
    #     print("Не удалось открыть файл", file=sys.stderr)
    #     return 1
    #
    # for i in range(N):
    #     file.write(f"{den[i]:.15g},")
    #
    # file.close()
    #
    # блок чтения, для возможности анализа весов скалярного произвдеения
    try:
        input = open(
            "/Users/bogdanvolkhin/Documents/Учеба ВУЗ/НГТУ/3 курс/ЧМ/ПЗ3/plot/den.txt",
            "r",
            encoding="utf-8",
        )
    except OSError:
        print("Не удалось открыть файл для чтения", file=sys.stderr)
        return 1

    den_read = []

    for val in input.read().split(","):
        if val.strip():
            den_read.append(float(val))

    input.close()

    # Для того чтобы получить обынчный интерполяционный сплайн нам достаточно указать параметр p = 0 и оставить вес w = 1, чтобы повторно не дублировать класс сплайна

    # Формируем вектор номеров точек, ля более удобной обработки
    x = [0.0] * N
    for i in range(N):
        x[i] = i + 1

    # cтроим интерполяцинный и сглажиывающий сплайны по одним и тем же данным

    interp = Smoothing_spline(0.0)

    smooth_1 = Smoothing_spline(0.4)

    smooth_2 = Smoothing_spline(0.8)

    smooth_3 = Smoothing_spline(0.99)

    interp.update_spl(x, den_read)
    smooth_1.update_spl(x, den_read)
    smooth_2.update_spl(x, den_read)
    smooth_3.update_spl(x, den_read)

    # формируем файл для сохранения значений сплайна в точках
    try:
        csv = open(
            "/Users/bogdanvolkhin/Documents/Учеба ВУЗ/НГТУ/3 курс/ЧМ/ПЗ3/plot/splines.csv",
            "w",
            encoding="utf-8",
        )
    except OSError:
        print("Не удалось открыть файл", file=sys.stderr)
        return 1

    # формируем соответсвующие значения в точках
    csv.write("x,interp,smooth_1,smooth_2,smooth_3\n")

    for i in range(N):
        csv.write(
            f"{x[i]:.15g},"
            f"{interp.get_spl_value_and_dev(x[i])[0]:.15g},"
            f"{smooth_1.get_spl_value_and_dev(x[i])[0]:.15g},"
            f"{smooth_2.get_spl_value_and_dev(x[i])[0]:.15g},"
            f"{smooth_3.get_spl_value_and_dev(x[i])[0]:.15g}\n"
        )
    csv.close()
    print("Данные записаны в splines.csv")

    return 0


if __name__ == "__main__":
    sys.exit(main(len(sys.argv), sys.argv))
