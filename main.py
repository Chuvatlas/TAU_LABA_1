import matplotlib.pyplot as pyplot
import control.matlab as matlab
import mpmath
import numpy as numpy
import sympy as sumpy
from math import *
import colorama as color


def graph(num, title, y, x):
    pyplot.subplot(2, 1, num)
    pyplot.tight_layout()  # отступы между графиками
    pyplot.grid(True)
    if title == 'АЧХ':
        pyplot.plot(x, y, 'green')
        pyplot.title(title)
        pyplot.ylabel('Magnitude')
        pyplot.xlabel('Omega (rad/s)')
    elif title == 'ФЧХ':
        pyplot.plot(x, y, 'yellow')
        pyplot.title(title)
        pyplot.ylabel('Phase (deg)')
        pyplot.xlabel('Omega (rad/s)')

def value_oscillation(h,t):
    t_p = matlab.stepinfo(w_closed, SettlingTimeThreshold=0.05)
    tp = t_p['SettlingTime']
    flg = False
    n = 0
    for h_i, t_j in zip(h, t):
        if t_j <= tp:
            delta = h_i- 1
            if delta > 0:
                if flg == False:
                    flg = True
                    n += 1
            if delta < 0:
                if flg == True:
                    flg = False
                    n += 1
    return n

def degree_of_attenuation(h, t):
    step_info = matlab.stepinfo(w_closed)
    flg = False
    n = 0
    h_m = []
    for h_i in h:
        if n < 3:
            delta = h_i - 1
            if delta > 0:
                if flg == False:
                    flg = True
                    n += 1
            if delta < 0:
                if flg == True:
                    flg = False
                    n += 1
        else:
            h_m.append(h_i)
            h2_max = max(h_m)
            psi = (step_info['Peak'] - h2_max)/(step_info['Peak'])
    return psi

def step_response(w_closed):
    print('Передаточная функция САУ : \n %s' % w_closed)
    pyplot.subplot()
    pyplot.grid(True)
    [h, t] = matlab.step(w_closed)
    pyplot.plot(t, h, 'red')
    pyplot.title('Переходная характеристика')
    pyplot.ylabel('Magnitude')
    pyplot.xlabel('Time (s)')
    pyplot.savefig('Step_response.png')
    pyplot.show()

    step_info = matlab.stepinfo(w_closed)
    print(step_info)
    t_p = matlab.stepinfo(w_closed, SettlingTimeThreshold=0.05)
    tp = t_p['SettlingTime']
    print('Время переходного процесса tп = ', tp)
    print('Перерегулирование: ', step_info['Overshoot'])
    print('Колебательность: ', value_oscillation(h, t))
    if value_oscillation(h, t) >= 3:
        print('Степень затухания: ', degree_of_attenuation(h, t))
    else:
        print('Степень затухания: 0')
    print('Величина первого максимума: ', step_info['Peak'], '\n'
          'Время достижения первого максимума: ', step_info['PeakTime'])
    #Интегральная оценка (численный метод)
    steady_state_value = step_info['SteadyStateValue']
    int_est = 0
    for i in range(1, len(t)):
        int_est += abs((h[i]-steady_state_value)*(t[i]-t[i-1]))
    print('Значение интеграла: ', int_est)

def roots_equation(w_closed):
    name = 'Корни Х.У.'
    print('Передаточная функция САУ : \n %s' % w_closed)
    pz = matlab.pzmap(w_closed)
    pyplot.grid(True)
    pyplot.show()
    pyplot.savefig('roots.png')

    poles = matlab.pole(w_closed)
    print("Полюса: \n %s" % poles)

    #Поиск времени регулирования
    re_max = numpy.max(numpy.real(poles))
    print('Время переходного процесса tп = ', abs(3/re_max))

    # Определение степени колебательности
    for comp_conj in poles:
        if (numpy.imag(comp_conj) > 0):
            alpha = numpy.real(comp_conj)
            beta = numpy.imag(comp_conj)
            print('beta = ', beta, '\n'
                  'alpha= ', alpha)
            mu = abs(beta/alpha)
    print('Степень колебалельности: ', mu)
    print('Перерегулирование: ', exp(-pi/mu))
    print('Степень затухания: ', 1-exp(-2*pi/mu))


def frequency_characteristics(w_closed):
    name = 'Диаграмма Боде'
    print('Передаточная функция разомкнутой САУ : \n %s' % w_opened)
    mag, phase, omega = matlab.bode(w_opened, dB=True)
    pyplot.show()
    gm, pm, wg, wp = matlab.margin(w_closed)
    print('Запас устойчиовости по амплитуде:', gm, ';\n',
          'Запас устойчивости по фазе:', pm, ';\n',
          'Критическая частота составляет:', wg, ';\n',
          'Частота среза составляет:', wp, '.\n')
    timeline = []
    for i in range(0, 1000):
        timeline.append(i / 1000)
    omega = []
    for i in range(0, 1000):
        omega.append(i / 1000)
    mag, phase, omega = matlab.freqresp(w_closed, timeline)
    graph(1, 'АЧХ', mag, omega)
    graph(2, 'ФЧХ', phase * 180 / pi, omega)

    print('Показатель колебательности: ', max(mag)/mag[0])
    print('Время переходного процесса: ', 1.5*2*pi/wp )
    pyplot.show()

# Начало кода здесь
type_regulator_select = input('Введите номер для выбора соответствующего регулятора : \n'
                  '1 - ПИ' + ';\n'
                  '2 - ПИД' + '.\n')
if type_regulator_select.isdigit():
    type_regulator_select = int(type_regulator_select)
    if type_regulator_select == 1:
        kp = 0.45 * (169 / 420)
        ki = 0.54 * ((169 / 420) / (2 * 20.9))
        w4 = matlab.tf([kp, ki], [1, 0])
    if type_regulator_select == 2:
        kp = 0.60 * (169 / 420)
        ki = 1.2 * ((169 / 420) / (2 * 20.9))
        kd = 0.075 * (169 / 420) * (2 * 20.9)
        w4 = matlab.tf([kd, kp, ki], [1, 0])
else:
    print(color.Fore.RED + '\nПожалуйста, введите числовое значение!')

w1 = matlab.tf([1], [8, 1])
w2 = matlab.tf([1], [5, 1])
w3 = matlab.tf([21], [5, 1])
# w4 = matlab.tf([169/420], [1])

w_closed = matlab.feedback(w1 * w2 * w3 * w4, 1)
print('Передаточная функция замкнутой САУ : \n %s' % w_closed)
w_opened = w1 * w2 * w3 * w4
print('Передаточная функция разомкнутой САУ : \n %s' % w_opened)

userInput = input('Введите номер характеристики : \n'
                  '1 - Переходная характеристика' + ';\n'
                  '2 - Корни Х.У.' + ';\n'
                  '3 - Часточные характериситики' + '.\n')
if userInput.isdigit():
    userInput = int(userInput)
    if userInput == 1:
        step_response = step_response(w_closed)
    elif userInput == 2:
        roots_equation = roots_equation(w_closed)
    elif userInput == 3:
        frequency_characteristics = frequency_characteristics(w_closed)
    else:
        print(color.Fore.RED + '\nНедопустимое числовое значение!')
else:
    print(color.Fore.RED + '\nПожалуйста, введите числовое значение!')
