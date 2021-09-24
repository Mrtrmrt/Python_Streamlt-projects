# -*- coding: utf-8 -*-
"""
Created on Wed Sep 22 15:50:13 2021

@author: Marat
"""

import numpy as np
import matplotlib.pyplot as plt
import streamlit as st
import math as m
import sympy as smp
from PIL import Image
import pandas as pd

#подрубаем интерфейс
st.title('Решение дифференциального уравнения конвекци-диффузии второй степени')
st.header('В данном файле представлена численная реализация решения уравнения конвекции-диффузии')
st.markdown('')
st.markdown('Исходное уравнение:')
st.latex(r'''
         \frac{\partial U}{\partial t}
         - \mu \cdot \frac{\partial{U^2}}{\partial{x^2}} = 0 
     ''')
st.markdown('Зададим пространственно-временную сетку:')
st.latex(r'''
         \widetilde \omega_{x} = \left\{x_{i} = ih, i = 0, 1, 2,...N \right\}
     ''')
st.latex(r'''
         \widetilde \omega_{\tau} = \left\{t_{j} = j\tau, j = 0, 1, 2,...M \right\}
     ''')
st.write(r"""
Подставим в соотвествие функции $U$ сеточную функцию  $y(i,j)$ и запишем разностные аналоги по координате и по времени:""")
st.latex(r'''
         \frac{\partial {U}}{\partial {t}}  = \frac{\partial {y}}{\partial {t}} 
         \sim \frac{y^{j+1}_{i} - y^{j}_{i}}{\tau}
        ''')
st.latex(r'''
         \frac{\partial {U}}{\partial {x}} = \frac{\partial {y}}{\partial {x}} 
         \sim \frac{y^{j}_{i+1} -2y^{j}_{i} + y^{j}_{i-1}}{h^2}
        ''')
st.write(r"""
Примем $\mu =1$, а начальное распределение зададим как функцию Гаусса:         
""")
st.latex(r'''
         U_{0}{(x)} = 0.75e^-\left({\frac{x-0.5}{0.1}} \right)^{2}
         ''')
st.markdown('Запишем итоговую разностную схему:')
st.latex(r'''
         \frac{y^{j+1}_{i} - y^{j}_{i}}{\tau} 
         - \mu \cdot \frac{y^{j}_{i+1} -2y^{j}_{i} + y^{j}_{i-1}}{h^2} = 0
     ''')
st.markdown('Итоговая разностная схема будет иметь вид:')

col1, col2, col3 = st.columns([3,6,1])
img = Image.open('E:\МФТИ\Вычислительная гидродинамика\ДЗ2\Scheme.png')
with col1:
    st.write("")

with col2:
    st.image(img)

with col3:
    st.write("")

st.markdown('Итоговое выражение для определения значения функции в точке на новом временном слое:')
st.latex(r'''
         y^{j+1}_{i} = \frac{\mu \tau}{h^2} \cdot
         (y^{j}_{i+1} -2y^{j}_{i} + y^{j}_{i-1}) + y^{j}_{i}
     ''')    
#Перейдем к решению
st.header('Решение')
st.subheader('Начальное распределение')


#блок констант
st.markdown('Константы:')
code_const = '''mu = 1'''
st.code(code_const)
mu = 1

#блок сетки по времени
st.markdown('Разбиение по времени:')
code_time = '''
tau = 0.001 #шаг по времени τ
t_start = 0 #старт
t_finish = 1 #финиш
t = np.arange(t_start,t_finish,tau)
'''
st.code(code_time)
tau = 0.001 #шаг по времени τ
t_start = 0
t_finish = 1
t = np.arange(t_start,t_finish + tau,tau)

#блок сетки по координате
st.markdown('Разбиение по координате:')
code_coord = '''
h = 0.05 #шаг по координате h
x_start = 0 #старт
x_finish = 1 #финиш
x = np.arange(x_start,x_finish+h,h)
'''
st.code(code_coord)
h = 0.05 #шаг по координате h
x_start = 0
x_finish = 1
x = np.arange(x_start,x_finish+h,h)

#блок исходной функции
st.markdown('Выведем график исходной функции:')
#определим значения исходной функции 
U = np.zeros(len(x))
for i in range(0,len(x)):
    #U.append(0.75*m.exp(-((x[i] - 0.5 )/0.5)**2))
    U[i] = 0.75*m.exp(-((x[i] - 0.5 )/0.1)**2)

#print(U)
print()

#функция для построения графика    
def Plot_U(coord, Z):
    
    fig, ax = plt.subplots(1,1, figsize = (10,8))
    ax.plot(coord, Z, 'c', alpha=1, linewidth=3)
    #Plot_U = plt.plot(coord, Z)
    plt.title('Начальная функция - Распределение Гаусса')
    plt.xlabel('x coordinate')
    plt.ylabel('U function')
    plt.grid()
    plt.show()
    st.pyplot(fig)
    
    return Plot_U

#вывод графика в streamlit
but_1 = st.checkbox('Вывести график')
#выделим место для кнопочек
if but_1:
    U_initial = Plot_U(x,U)
    
U0 = U

#функция решения для аппроксимации второй производной через центральную разность
def Solve_C(time_C, coord_C, U_init_C):
    
    C_C = len(coord_C)
    T_C = len(time_C)
    
    dt_C = round(time_C[-1] / T_C, 3) #вычисляем шаг по времени для переданного массива
    dx_C = round(coord_C[-1] / C_C, 3) #вычисляем шаг по координате для переданного массива
    buf = ((mu * dt_C)/ (dx_C)**2)
    
    Sol_new_C = np.zeros(C_C) #массив итогового решения
    Sol_old_C = np.zeros(C_C) #состояние системы на новом временном слое
    
    for k in range(0,len(U_init_C)):
        Sol_new_C[k] = U_init_C[k]
 
    #j - цикл по времени, поскольку ищем решение на новом временном слое
    for j in range(0, T_C):
        Sol_old_C = Sol_new_C.copy()#после завершения каждой итерации в этот массив будут сохраняться
        #print(Sol_new_C)
        #i - цикл по координате x
        for i in range(1, C_C-1):
                Sol_new_C[i] = buf * (Sol_old_C[i+1] - 2 * Sol_old_C[i] + Sol_old_C[i-1]) + Sol_old_C[i]

    return Sol_new_C

#функция построения графика полученного решения
def Plotting(coord, Solution, Name, time_step):
    
    fig, ax = plt.subplots(1,1, figsize = (10,8))
    ax.plot(coord, Solution, 'c', alpha=1, linewidth=3)
    plt.title('Полученное решение при t = {} c, методом {}'.format(time_step, Name))
    plt.xlabel('x coordinate')
    plt.ylabel('Решение')
    plt.grid()
    plt.show()
    st.pyplot(fig)
    
    return Plotting

#блок полученных решений

#решение в нулевой момент времени
st.subheader('Решение в момент времени t =0')
st.markdown('Полученный график в данном случае будет аналогичен исходной функции:')

#выделим место для кнопочек
but_2 = st.checkbox('Вывести график для t=0')
#условие нажатия
if but_2:
    U_initial = Plot_U(x,U)

#решение в момент времени t=0.05
st.subheader('Решение в момент времени t=0.05 c')
st.markdown('Выведем график решения задачи в момент времени $t = 0.05 c$ :')

#Вызов расчетных функций
Method = 'FTCS'
t_2 = 0.05
i_t_2 = int(t_2 / tau)
#теперь вызовем функцию Solve
Sol_t_2 = Solve_C(t[0:i_t_2], x, U0)

#выделим место для кнопочек
but_3 = st.checkbox('Вывести график для t=0.05 c')
#условие нажатия
if but_3:
    Plot_second = Plotting(x, Sol_t_2, Method, t_2)

#решение в момент времени t=0.01
st.subheader('Решение в момент времени t=0.1 c')
st.markdown('Выведем график решения задачи в момент времени $t = 0.1 c$ :')

#Вызов расчетных функций
t_3 = 0.1
i_t_3 = int(t_3 / tau)
#теперь вызовем функцию Solve
Sol_t_3 = Solve_C(t[0:i_t_3], x, U0)

#выделим место для кнопочек
but_4 = st.checkbox('Вывести график для t=0.1 c')
#условие нажатия
if but_4:
    Plot_third = Plotting(x, Sol_t_3, Method, t_3)