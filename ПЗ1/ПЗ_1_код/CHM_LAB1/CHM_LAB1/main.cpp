#include <cmath>
#include <iomanip>
#include <iostream>
#include <stdexcept>
#include <vector>


using namespace std;

// Граница для сравнения с еденицей
double eps_sr = 1e-14;
// Граница точности по условию задачи
double eps = 0.0001;


//-------------------------------------------------------- ПУНКТ 1
// Опишем функции для регулярного и адаптивного сеточного разбиения

// Функция регулярного разбиения сетки
vector<double> generateRegGrid(const double &a, const double &b, int n) {
    // Проверка на факт того, что количество участков разбиений отрезка положительно
    if (n <= 0) {
        // Генерируем исключение
        throw invalid_argument("Количество сегментов должно быть положительным.");
    }

    // Определяем динамический массив длиной n+1
    vector<double> x(n + 1);

    // Определяем размер одного сегмента
    double h = (b - a) / n;

    // Формируем узлы сетки
    for (int i = 0; i <= n; ++i) {
        if(i != n) {
            x[i] = a + i * h;
        }else {
            x[i] = b;
        }

    }

    return x;
}

// Функция адаптивного разбиения сетки
vector<double> generateAdaptiveGrid(const double &a, const double &b, int n, double r) {
    // Проверка на факт того, что количество участков разбиений отрезка положительно
    if (n <= 0) {
        throw invalid_argument("Количество сегментов должно быть положительным.");
    }

    // Проверка на факт того, что r число >=0
    if (r <= 0) {
        throw invalid_argument("Коэффициент r должен быть больше нуля.");
    }

    // Определяем динамический массив длиной n+1
    vector<double> x(n + 1);

    // Задааем первую точку разбиения как начало отрезка
    x[0] = a;

    // Если r = 1, получаем обычную равномерную сетку
    if (abs(r - 1.0) < eps_sr) {
        return generateRegGrid(a, b, n);
    }

    // Для поиска h1 используем формлу геометрической прогресси
    double h1 = (b - a) * (1.0 - r) /
                (1.0 - pow(r, n));

    // Текущий шаг
    double h = h1;

    // Построение узлов
    for (int i = 1; i <= n; ++i) {
        x[i] = x[i - 1] + h;

        // Следующий шаг отличается от предыдущего в r раз
        h *= r;
    }

    x[n] = b;

    return x;
}

// Функция печати сетки
void printGrid(const vector<double> &x) {
    //Фиксируем вывод чисел с плавающей точкой с 6 занками после запятой
    cout << fixed << setprecision(6);

    //Выводим заголовок таблицы
    cout << "\n";
    cout << "---------------------------------------------\n";
    cout << " i\t\tx_i\t\t h_i\n";
    cout << "---------------------------------------------\n";

    // В цикле выводим соответсвующие значения x и h построчно
    for (size_t i = 0; i < x.size(); ++i) {
        cout << i << "\t\t" << x[i];

        // Для последнего узла следующего шага нет
        if (i + 1 < x.size()) {
            double h = x[i + 1] - x[i];
            cout << "\t\t" << h;
        }

        cout << "\n";
    }

    cout << "---------------------------------------------\n";
}

//-------------------------------------------------------- ПУНКТ 2

// Класс реализующий интерфейс кубического интерполяционного сплайна
class CubicSpline {
    private:
        //Набор точек и значений табличной функции для формирования сплайна
        vector<double> x, y;

        // Коэффициенты сплайна
        vector<double> a, b, c, d;

        // Длины сегментов
        vector<double> h;

        // Коэффициенты трехдиагональной системы
        vector<double> alpha, beta, gamma, phi;

        // Коэффициенты прямого хода
        vector<double> sigma, psi;

    public:

        // Констурктор класса
        CubicSpline(const vector<double> &x_val, const vector<double> &y_val) {
            build(x_val, y_val);
        }

        // Функция формирования кубического спалйна
        void build(const vector<double> &x_val, const vector<double> &y_val) {

            // Проверяем, что количество точек по x и по y сопадают
            if (x_val.size() != y_val.size()) {
                throw invalid_argument("Количество x и y должно совпадать.");
            }

            // Проверяем, что количество точек миниум две, чтобы был хотя-бы один сегмент
            if (x_val.size() < 2) {
                throw invalid_argument("Необходимо минимум 2 узла.");
            }

            // Определяем количество сегментов, которое равно числу точек -1
            int n = static_cast<int>(x_val.size()) - 1;

            // Сохраняем исходные данные
            x = x_val;
            y = y_val;

            // Выделяем память под соответсвтующие векторы
            a.resize(n + 1);
            b.resize(n + 1);
            c.resize(n + 1);
            d.resize(n + 1);

            h.resize(n);

            /*
               alpha, beta, gamma, phi,
               sigma, psi используются только
               для внутренних узлов:

               Поэтому достаточно n элементов.
             */

            alpha.resize(n);
            beta.resize(n);
            gamma.resize(n);
            phi.resize(n);

            sigma.resize(n);
            psi.resize(n);

            // Вычисляем длины всех сегментов
            for (int k = 0; k < n; ++k) {
                h[k] = x[k + 1] - x[k];
            }

            // Формируем трехдиагональную систему
            for (int k = 1; k <= n - 1; ++k) {

                alpha[k] = h[k - 1];

                beta[k] = 2.0 * (h[k - 1] + h[k]);

                gamma[k] = h[k];

                // Правая часть phi_k
                phi[k] = 3.0 * ((y[k + 1] - y[k]) / h[k] - (y[k] - y[k - 1]) / h[k - 1]);
            }

            // Используем граничные условия
            c[0] = 0.0;
            c[n] = 0.0;

            // Осуществляем прямой ход
            if (n >= 2) {
                // Элемнты первой строки не изменяются
                sigma[1] = beta[1];
                psi[1] = phi[1];

                for (int k = 2; k <= n - 1; ++k) {
                    sigma[k] = beta[k] - alpha[k] * gamma[k - 1] / sigma[k - 1];

                    psi[k] = phi[k] - alpha[k] * psi[k - 1] / sigma[k - 1];
                }

                // Обратный ход

                c[n - 1] = psi[n - 1] / sigma[n - 1];

                for (int k = n - 2; k >= 1; --k) {
                    c[k] = (psi[k] - gamma[k] * c[k + 1]) / sigma[k];
                }
            }

            for (int k = 0; k < n; ++k) {
                a[k] = y[k];

                d[k] = (c[k + 1] - c[k]) / (3.0 * h[k]);

                b[k] = (y[k + 1] - y[k]) / h[k] - h[k] * (2.0 * c[k] + c[k + 1]) / 3.0;
            }
        }

        // Функция вычисления S(x)
        double calculate(double X) const {
            // Формируем число сегментов
            int n = static_cast<int>(x.size()) - 1;

            // Проверяем диапазон
            if (X < x[0] || X > x[n]) {
                throw out_of_range("Точка X находится за пределами отрезка.");
            }

            // Определяем необходимый сегмент

            int k = 0;

            for (int i = 0; i < n; ++i) {
                if (X >= x[i] && X <= x[i + 1]) {
                    k = i;
                    break;
                }
            }

            // Вычисляем локальную координату
            double t = X - x[k];

            // Вычисляем значение кубического полинома
            return a[k] + b[k] * t + c[k] * t * t + d[k] * t * t * t;
        }

        // Вычисление S'(x)
        double calculate_first_der(double X) const {

            // Формируем число сегментов
            int n = static_cast<int>(x.size()) - 1;

            // Проверяем диапазон
            if (X < x[0] || X > x[n]) {
                throw out_of_range("Точка X находится за пределами отрезка.");
            }

            // Определяем необходимый сегмент

            int k = 0;

            for (int i = 0; i < n; ++i) {
                if (X >= x[i] && X <= x[i + 1]) {
                    k = i;
                    break;
                }
            }

            // Вычисляем локальную координату
            double t = X - x[k];

            // Вычисляем значение кубического полинома
            return b[k] + 2 * c[k] * t + 3 * d[k] * t * t;
        }

        // Вычисление S''(x)
        double calculate_second_der(double X) const {

            // Формируем число сегментов
            int n = static_cast<int>(x.size()) - 1;

            // Проверяем диапазон
            if (X < x[0] || X > x[n]) {
                throw out_of_range("Точка X находится за пределами отрезка.");
            }

            // Определяем необходимый сегмент

            int k = 0;

            for (int i = 0; i < n; ++i) {
                if (X >= x[i] && X <= x[i + 1]) {
                    k = i;
                    break;
                }
            }

            // Вычисляем локальную координату
            double t = X - x[k];

            // Вычисляем значение кубического полинома
            return 2 * c[k] + 6 * d[k] * t;
        }

        // Функция вычисления точности аппроксимации
        double calculate_error(const double &a, const double &b, int derivative) const {
            // Проходимся по всем элемнтам массива точек и вычисляем значение функции в этих точкаx

            const int N = 100000;
            double max_error = 0.0;

            double h = (b - a) / N;

            for (int i = 0; i <= N; i++) {
                double fun;
                double spl;
                double X;
                if (i == N)
                    X = b;
                else
                    X = a + h * i;


                if (derivative == 0) {
                    fun = sin(X);
                    spl = calculate(X);
                } else if (derivative == 1) {
                    fun = cos(X);
                    spl = calculate_first_der(X);
                } else if (derivative == 2) {
                    fun = -sin(X);
                    spl = calculate_second_der(X);
                } else {
                    throw invalid_argument("Неверный порядок производной.");
                }

                // Вычисляем модуль разности значений сплайна и функции
                double error = abs(spl - fun);

                // Вычисляем максимальную ошибку
                if (error > max_error) {
                    max_error = error;
                }
            }

            return max_error;
        }

        // Функция вывода коэффициентов сплайна
        void printCoefficients() const {
            
            int n = static_cast<int>(x.size()) - 1;

            cout << fixed << setprecision(6);

            cout << "\n";
            cout << "====================================================\n";
            cout << "        КОЭФФИЦИЕНТЫ КУБИЧЕСКОГО СПЛАЙНА\n";
            cout << "====================================================\n";

            cout << "\n";

            cout << " k"
                 << "\t a_k"
                 << "\t\t b_k"
                 << "\t\t c_k"
                 << "\t\t d_k\n";

            cout << "----------------------------------------------------\n";

            for (int k = 0; k < n; ++k) {
                cout << k
                     << "\t " << a[k]
                     << "\t " << b[k]
                     << "\t " << c[k]
                     << "\t " << d[k]
                     << "\n";
            }

            cout << "----------------------------------------------------\n";
        }

        // Функция вывода коэффициентов трехдиагональной системы
        void printSystemCoefficients() const {
            int n = static_cast<int>(x.size()) - 1;

            cout << fixed << setprecision(6);

            cout << "\n";
            cout << "Коэффициенты системы:\n\n";

            cout << " k"
                 << "\t alpha_k"
                 << "\t beta_k"
                 << "\t gamma_k"
                 << "\t phi_k\n";

            cout << "----------------------------------------------------\n";

            for (int k = 1; k <= n - 1; ++k) {
                cout << k
                     << "\t " << alpha[k]
                     << "\t " << beta[k]
                     << "\t " << gamma[k]
                     << "\t " << phi[k]
                     << "\n";
            }

            cout << "\n";
            cout << "Коэффициенты после прямого хода:\n\n";

            cout << " k"
                 << "\t sigma_k"
                 << "\t psi_k\n";

            cout << "----------------------------------------------------\n";

            for (int k = 1; k <= n - 1; ++k) {
                cout << k
                     << "\t " << sigma[k]
                     << "\t " << psi[k]
                     << "\n";
            }

            cout << "----------------------------------------------------\n";
        }

};


//-------------------------------------------------------- ПУНКТ 3

int main() {
    double a, b;
    int n;

    // Осуществляем ввод границ отрезка и числа элемнтов разбиения
    cout << "Введите a и b: ";
    cin >> a >> b;
    cout << "Введите количество сегментов, на которые необходимо осуществлять равномерное разбиение: ";
    cin >> n;

    // Вводим дополнительное условие на число сегментов, чтобы выполнялось условие на число точек из задания
    if (n < 10) {
        throw invalid_argument("Количество сегментов должно быть не менее 10.");
    }
    
    // формируем равномерное сеточное разбиение
    vector<double>x = generateRegGrid(a, b, n);
    vector<double>y(x.size());

    // в качестве непрерывной неполиномиальной функции рассмотрим функцию синуса (на вход принимаются радианы, а не градусы)
    // формируем табличную функцию по значениям в узлах сетки
    for (int i = 0; i < x.size(); i++) {
        y[i] = sin(x[i]);
    }

    // Сформируем таблицу значений сплайна и его двух первых производных в точках, которые не совпадают с узловыми

    // Создаем кубический сплайн
    CubicSpline spline(x, y);
    // Выводим его коэффициенты
    spline.printCoefficients();

    double ln = (b - a) / n;

    vector<double>x_not_base(n);
    vector<double>spline_out(n);
    vector<double>spline_out_dev(n);
    vector<double>spline_out_dev2(n);

    for (int i = 0; i < n; i++) {
        double k = 0.1 + 0.8 * ((double)rand() / RAND_MAX);

        x_not_base[i] = a + k * (ln) + ln * i;

        spline_out[i] = spline.calculate(x_not_base[i]);
        spline_out_dev[i] = spline.calculate_first_der(x_not_base[i]);
        spline_out_dev2[i] = spline.calculate_second_der(x_not_base[i]);
    }

    cout << "\nТаблица значений сплайна:\n";
    cout << "x\t\tS(x)\t\tS'(x)\t\tS''(x)\n";

    for (int i = 0; i < n; ++i) {
        cout << fixed << setprecision(8)
             << x_not_base[i] << "\t"
             << spline_out[i] << "\t"
             << spline_out_dev[i] << "\t"
             << spline_out_dev2[i] << "\n";
    }

    double error = spline.calculate_error(a, b, 0);
    double error_dev = spline.calculate_error(a, b, 1);
    double error_dev2 = spline.calculate_error(a, b, 2);

    cout << "Погрешность S(x):   " << error << endl;
    cout << "Погрешность S'(x):  " << error_dev << endl;
    cout << "Погрешность S''(x): " << error_dev2 << endl;
    //----------------------------------------------
    vector<double>x_h_2 = generateRegGrid(a, b, n * 2);

    vector<double>y_h_2(x_h_2.size());

    for (int i = 0; i < x_h_2.size(); i++) {
        y_h_2[i] = sin(x_h_2[i]);
    }

    // Создаем кубический сплайн
    CubicSpline spline_h_2(x_h_2, y_h_2);
    // Выводим его коэффициенты
    spline_h_2.printCoefficients();

    double ln_h_2 = (b - a) / (n * 2);

    vector<double>x_not_base_h_2(2 * n);
    vector<double>spline_out_h_2(2 * n);
    vector<double>spline_out_dev_h_2(2 * n);
    vector<double>spline_out_dev2_h_2(2 * n);

    for (int i = 0; i < 2 * n; i++) {
        double k = 0.1 + 0.8 * ((double)rand() / RAND_MAX);

        x_not_base_h_2[i] = a + k * (ln_h_2) + ln_h_2 * i;

        spline_out_h_2[i] = spline_h_2.calculate(x_not_base_h_2[i]);
        spline_out_dev_h_2[i] = spline_h_2.calculate_first_der(x_not_base_h_2[i]);
        spline_out_dev2_h_2[i] = spline_h_2.calculate_second_der(x_not_base_h_2[i]);
    }

    cout << "\nТаблица значений сплайна:\n";
    cout << "x\t\tS(x)\t\tS'(x)\t\tS''(x)\n";

    for (int i = 0; i < 2 * n; ++i) {
        cout << fixed << setprecision(8)
             << x_not_base_h_2[i] << "\t"
             << spline_out_h_2[i] << "\t"
             << spline_out_dev_h_2[i] << "\t"
             << spline_out_dev2_h_2[i] << "\n";
    }

    double error_h_2 = spline_h_2.calculate_error(a, b, 0);
    double error_dev_h_2 = spline_h_2.calculate_error(a, b, 1);
    double error_dev2_h_2 = spline_h_2.calculate_error(a, b, 2);

    cout << "Погрешность S(x):   " << error_h_2 << endl;
    cout << "Погрешность S'(x):  " << error_dev_h_2 << endl;
    cout << "Погрешность S''(x): " << error_dev2_h_2 << endl;
    //----------------------------------------------

    vector<double>x_h_4 = generateRegGrid(a, b, n * 4);

    vector<double>y_h_4(x_h_4.size());

    for (int i = 0; i < x_h_4.size(); i++) {
        y_h_4[i] = sin(x_h_4[i]);
    }

    // Создаем кубический сплайн
    CubicSpline spline_h_4(x_h_4, y_h_4);
    // Выводим его коэффициенты
    spline_h_4.printCoefficients();

    double ln_h_4 = (b - a) / (n * 4);

    vector<double>x_not_base_h_4(4 * n);
    vector<double>spline_out_h_4(4 * n);
    vector<double>spline_out_dev_h_4(4 * n);
    vector<double>spline_out_dev2_h_4(4 * n);

    for (int i = 0; i < 4 * n; i++) {
        double k = 0.1 + 0.8 * ((double)rand() / RAND_MAX);

        x_not_base_h_4[i] = a + k * (ln_h_4) + ln_h_4 * i;

        spline_out_h_4[i] = spline_h_4.calculate(x_not_base_h_4[i]);
        spline_out_dev_h_4[i] = spline_h_4.calculate_first_der(x_not_base_h_4[i]);
        spline_out_dev2_h_4[i] = spline_h_4.calculate_second_der(x_not_base_h_4[i]);
    }

    cout << "\nТаблица значений сплайна:\n";
    cout << "x\t\tS(x)\t\tS'(x)\t\tS''(x)\n";

    for (int i = 0; i < 4 * n; ++i) {
        cout << fixed << setprecision(8)
             << x_not_base_h_4[i] << "\t"
             << spline_out_h_4[i] << "\t"
             << spline_out_dev_h_4[i] << "\t"
             << spline_out_dev2_h_4[i] << "\n";
    }

    double error_h_4 = spline_h_4.calculate_error(a, b, 0);
    double error_dev_h_4 = spline_h_4.calculate_error(a, b, 1);
    double error_dev2_h_4 = spline_h_4.calculate_error(a, b, 2);

    cout << "Погрешность S(x):   " << error_h_4 << endl;
    cout << "Погрешность S'(x):  " << error_dev_h_4 << endl;
    cout << "Погрешность S''(x): " << error_dev2_h_4 << endl;
    
    
    //-------------------------------------------------------- ПУНКТ 4

    // Блок сравнения сплайна и метода конечных разностей

    vector<double> first_disc_dev_spl(x.size());
    vector<double> second_disc_dev_spl(x.size());

    // Вычисляем значения  первой и второй производной в точках исходного равноменрного разбиения при использовании сплайна
    for (int i = 0; i < x.size(); i++) {
        first_disc_dev_spl[i] = spline.calculate_first_der(x[i]);
        second_disc_dev_spl[i] = spline.calculate_second_der(x[i]);
    }

    double max_error_frst_spl = 0.0;
    for (int i = 0; i < x.size(); i++) {

        double fun = cos(x[i]);

        // Вычисляем ошибку на конкретной итерации
        double error = abs(first_disc_dev_spl[i] - fun);

        // Вычисляем максимальную ошибку
        if (error > max_error_frst_spl) {
            max_error_frst_spl = error;
        }

    }

    double max_error_sec_spl = 0.0;
    for (int i = 0; i < x.size(); i++) {

        double fun = -sin(x[i]);

        // Вычисляем ошибку на конкретной итерации
        double error = abs(second_disc_dev_spl[i] - fun);

        // Вычисляем максимальную ошибку
        if (error > max_error_sec_spl) {
            max_error_sec_spl = error;
        }

    }

    cout << "Погрешность S'(x) disc:  " << max_error_frst_spl << endl;
    cout << "Погрешность S''(x) disc: " << max_error_sec_spl << endl;


    // В рамках вычисления дискретной производной функции методом коненчных разностей будем рассматривать метод численного дифференцирования 2-го порядка точности.
    // Поскольку мы используем в рассчетах равномерное разбиение, то формулы производных 1-го и 2-го порядка примут вид

    vector<double> first_disc_dev(x.size());
    vector<double> second_disc_dev(x.size());

    // Вычислеяем значения во внутренних узлах
    for (int i = 1; i < x.size() - 1; i++) {
        first_disc_dev[i] = (sin(x[i + 1]) - sin(x[i - 1])) / (2 * (x[1] - x[0]));
        second_disc_dev[i] = ((sin(x[i + 1]) - 2 * sin(x[i]) + sin(x[i - 1])) / ((x[1] - x[0]) * (x[1] - x[0])));
    }
    // Вычисляем значения в краевых узлах
    first_disc_dev[0] = (-3 * sin(x[0]) + 4 * sin(x[1]) - sin(x[2])) / (2 * (x[1] - x[0]));
    first_disc_dev[x.size() - 1] = (3 * sin(x[x.size() - 1]) - 4 * sin(x[x.size() - 2]) + sin(x[x.size() - 3])) /
                                   (2 * (x[1] - x[0]));

    second_disc_dev[0] = (2 * sin(x[0]) - 5 * sin(x[1]) + 4 * sin(x[2]) - sin(x[3])) / ((x[1] - x[0]) * (x[1] - x[0]));
    second_disc_dev[x.size() - 1] = (2 * sin(x[x.size() - 1]) - 5 * sin(x[x.size() - 2]) + 4 * sin(x[x.size() - 3]) -
                                     sin(x[x.size() - 4])) / ((x[1] - x[0]) * (x[1] - x[0]));

    double max_error_frst = 0.0;
    for (int i = 0; i < x.size(); i++) {

        double fun = cos(x[i]);

        // Вычисляем ошибку на конкретной итерации
        double error = abs(first_disc_dev[i] - fun);

        // Вычисляем максимальную ошибку
        if (error > max_error_frst) {
            max_error_frst = error;
        }

    }

    double max_error_sec = 0.0;
    for (int i = 0; i < x.size(); i++) {

        double fun = -sin(x[i]);

        // Вычисляем ошибку на конкретной итерации
        double error = abs(second_disc_dev[i] - fun);

        // Вычисляем максимальную ошибку
        if (error > max_error_sec) {
            max_error_sec = error;
        }

    }

    cout << "Погрешность f'(x) disc:  " << max_error_frst << endl;
    cout << "Погрешность f''(x) disc: " << max_error_sec << endl;


    return 0;
}
