#include <cmath>
#include <iomanip>
#include <iostream>
#include <stdexcept>
#include <vector>
#include <functional>
#include <numbers>

#include <random>

#include <fstream>

using namespace std;

// Опишем класс реализующий интерфейс сглаживающего сплайна

class Smoothing_spline {

    private:
        // параметр сглаживания
        double SMOOTH;
        // точки сетки
        vector<double> x;
        // коэффициенты разложения по базису
        vector<double> alpha;
        //
        vector<double> a, b, c;

        // перход на мастер элемнент [-1,1]
        double Trans_To_Mas_El(int s_num, const double &X){
            double Ksi = 2.0 * (X - x[s_num]) / (x[s_num + 1] - x[s_num]) - 1.0;
            return Ksi;
        }

        // базисные функции на [-1,1]
        double Bas_Func(int N, const double &Ksi){
            switch (N) {
                case 1:
                {
                    return 0.5 * (1 - Ksi);
                    break;
                }

                case 2:
                {
                    return 0.5 * (1 + Ksi);
                    break;
                }

                default: {
                    throw invalid_argument("Укажите корректный номер функции");
                    break;
                }

            }
        }

        // производные базисных функций на [-1,1]
        double Der_Bas_Func(int N, const double &Ksi){

            switch (N) {
                case 1:
                {
                    return -0.5;
                    break;
                }

                case 2:
                {
                    return 0.5;
                    break;
                }

                default: {
                    throw invalid_argument("Укажите корректный номер функции");
                    break;
                }

            }

        }

    public:
        // конструктор сглаживающего сплайна
        Smoothing_spline(const double &SMOOTH){
            this->SMOOTH = SMOOTH;
        }

        // Функция для обновления сплайна
        void update_spl(const vector<double> &x_upd, const vector<double> &f_d){

            // обновляем точки сетки
            x.clear();

            for(auto &elem :x_upd) x.push_back(elem);

            // определяем число отрезков разбиения
            int n = static_cast<int>(x.size()) - 1;

            // формируем вектор коэффициентов разложения по базису
            alpha.resize(n + 1);

            // описываем диагонали матрицы СЛАУ

            a.resize(n + 1);
            b.resize(n + 1);
            c.resize(n + 1);

            // процедура ассемблирования СЛАУ
            auto Assembly = [&](int i, double X, double F_Val, double w) {

                                // переход на мастер элемент
                                double Ksi = Trans_To_Mas_El(i, X);
                                // вычисление значений баззисных функций на мастер элементе
                                double f1 = Bas_Func(1, Ksi);
                                double f2 = Bas_Func(2, Ksi);
                                // внесение вкладов в матрицу СЛАУ
                                b[i] += (1.0 - SMOOTH) * w * f1 * f1;
                                b[i + 1] += (1.0 - SMOOTH) * w * f2 * f2;

                                a[i + 1] += (1.0 - SMOOTH) * w * f1 * f2;
                                c[i] += (1.0 - SMOOTH) * w * f2 * f1;
                                // внесение вкладов в вектор правой части СЛАУ
                                alpha[i] += (1.0 - SMOOTH) * w * f1 * F_Val;
                                alpha[i + 1] += (1.0 - SMOOTH) * w * f2 * F_Val;
                            };
            // сборка СЛАУ
            for (int i = 0; i < n; i++) {

                double W_left  = (i == 8)     ? 0.5 : 1.0;
                double W_right = (i + 1 == 8) ? 0.5 : 1.0;

                Assembly(i, this->x[i], f_d[i], W_left);
                Assembly(i, this->x[i + 1], f_d[i + 1], W_right);

                double h = x[i + 1] - x[i];
                b[i] += 1.0 / h * SMOOTH;
                b[i + 1] += 1.0 / h * SMOOTH;
                a[i + 1] -= 1.0 / h * SMOOTH;
                c[i] -= 1.0 / h * SMOOTH;

            }

            // метод прогонки: прямой ход
            for(int j = 1; j < n + 1; j++) {
                b[j] -= a[j] / b[j - 1] * c[j - 1];
                alpha[j] -= a[j] / b[j - 1] * alpha[j - 1];
            }
            // метод прогонки: обратный ход
            alpha[n] /= b[n];
            for(int j = n - 1; j >= 0; j--) {
                alpha[j] = (alpha[j] - alpha[j + 1] * c[j]) / b[j];
            }


        }

        // Функция для вычисления значения сплайна и его производной в точке
        vector<double> get_spl_value_and_dev(const double &x_in){

            // машинный ноль
            double eps = 1e-7;

            // определяем число отрезков разбиения
            int n = static_cast<int>(x.size()) - 1;

            // вектор результата
            vector<double> res(3);

            for(int i = 0; i < n; i++) {
                if((x_in > x[i] && x_in < x[i + 1]) || fabs(x_in - x[i]) < eps || fabs(x_in - x[i + 1]) < eps) {
                    // длина отрезка
                    double local = x[i + 1] - x[i];
                    // переход на мастер элемент
                    double Ksi;
                    Ksi = Trans_To_Mas_El(i, x_in);
                    // вычисление значения сплайна и его производной по базисным функциям
                    res[0] = alpha[i] * Bas_Func(1, Ksi) + alpha[i + 1] * Bas_Func(2, Ksi);
                    res[1] = (alpha[i] * Der_Bas_Func(1, Ksi) + alpha[i + 1] * Der_Bas_Func(2, Ksi)) * 2.0 / local;
                    res[2] = 0.0;

                    return res;
                }

            }
            throw invalid_argument("Указано не корректное число");
        }

};


// Генерация случайных чисел из нормального распределения
vector<double> normDen(double M, double s, const int N)
{
    if (!isfinite(M)) {
        throw invalid_argument(
                  "Мат ожидание должно быть конечным числом");
    }
    if (!isfinite(s) || s <= 0.0) {
        throw invalid_argument(
                  "Стандартное отклонение должно быть конечным и положительным");
    }

    static mt19937 gen(random_device{}());
    normal_distribution<double> dist(M, s);

    vector<double> val;
    val.reserve(N);

    for (int i = 0; i < N; ++i) {
        val.push_back(dist(gen));
    }

    return val;
}

int main(int argc, const char *argv[]) {

    // параметры нормального распредления и число точек

    int N = 1443;
    double M = 0.94;
    double s = 4.28;

    // формируем последовательность случайных чисел по нормальному закону

    // блок записи, для возможности проведения анализа весов скалярного произвдеения
    //const vector<double>den = normDen(M, s, N);

//    ofstream file(
//        "/Users/bogdanvolkhin/Documents/Учеба ВУЗ/НГТУ/3 курс/ЧМ/ПЗ3/plot/den.txt"
//    );
//    if (!file.is_open()) {
//        cerr << "Не удалось открыть файл\n";
//        return 1;
//    }
//
//    file << setprecision(15);
//
//    for (int i = 0; i < N; ++i) {
//        file << den[i] << ',';
//    }
//
//    file.close();
//
    // блок чтения, для возможности анализа весов скалярного произвдеения
    ifstream input(
        "/Users/bogdanvolkhin/Documents/Учеба ВУЗ/НГТУ/3 курс/ЧМ/ПЗ3/plot/den.txt"
        );

    if (!input.is_open()) {
        cerr << "Не удалось открыть файл для чтения\n";
        return 1;
    }

    vector<double> den_read;
    double val;
    char sep;

    while (input >> val) {
        den_read.push_back(val);
        if (!(input >> sep)) {
            break;
        }
        if (sep != ',') {
            cerr << "Ожидалась запятая\n";
            return 1;
        }
    }

    input.close();

    // Для того чтобы получить обынчный интерполяционный сплайн нам достаточно указать параметр p = 0 и оставить вес w = 1, чтобы повторно не дублировать класс сплайна

    // Формируем вектор номеров точек, для более удобной обработки
    vector<double> x(N);
    for (int i = 0; i < N; ++i) {
        x[i] = i + 1;
    }

    // cтроим интерполяцинный и сглажиывающий сплайны по одним и тем же данным

    Smoothing_spline interp(0.0);

    Smoothing_spline smooth_1(0.4);

    Smoothing_spline smooth_2(0.8);

    Smoothing_spline smooth_3(0.99);

    interp.update_spl(x, den_read);
    smooth_1.update_spl(x, den_read);
    smooth_2.update_spl(x, den_read);
    smooth_3.update_spl(x, den_read);

    // формируем файл для сохранения значений сплайна в точках
    ofstream csv(
        "/Users/bogdanvolkhin/Documents/Учеба ВУЗ/НГТУ/3 курс/ЧМ/ПЗ3/plot/splines.csv"
        );

    if (!csv.is_open()) {
        cerr << "Не удалось открыть файл\n";
        return 1;
    }

    csv << setprecision(15);
    // формируем соответсвующие значения в точках
    csv << "x,interp,smooth_1,smooth_2,smooth_3\n";

    for (int i = 0; i < N; ++i) {
        csv << x[i] << ','
            << interp.get_spl_value_and_dev(x[i])[0] << ','
            << smooth_1.get_spl_value_and_dev(x[i])[0] << ','
            << smooth_2.get_spl_value_and_dev(x[i])[0] << ','
            << smooth_3.get_spl_value_and_dev(x[i])[0] << '\n';
    }
    csv.close();
    cout << "Данные записаны в splines.csv\n";

    return 0;
}
