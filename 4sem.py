import matplotlib.pyplot as plt
import numpy as np

from numpy import log, sqrt, exp
from numpy import cos, sin, pi
from random import random
from mpl_toolkits.mplot3d import Axes3D

def square(x):
    return x * x

# количество столкновений
n = 1000

# количество электронов (K)
k = 200

# лямбда
lbd = 1

# скорость электрона
v = 1

# двумерный массив времен столкновений
t_array = np.zeros((3, n))

# трёхмерный массив с квадратами смещений (по каждому из K электронов для n столкновений)
l_square_array = np.zeros((3, k, n))

# двумерный массив усреднённых квадратов смещений
average_l_square_array = np.zeros((3, n))

# массивы с координатами K электронов в n столкновениях (3 случая)
x_array = np.zeros((3, k, n))
y_array = np.zeros((3, k, n))
z_array = np.zeros((3, k, n))

# перебор всех электронов
for i in range(k):
    
    cos_t_3 = np.zeros(3)
    r = random()
    cos_t_3[0] = 1 - 2 * r
    r = random()
    cos_t_3[1] = 1 - 2 * pow(r, 2/3)
    r = random()
    cos_t_3[2] = 2 * pow(1 - r, 2/3) - 1
    
    l_3 = np.zeros(3)
    fi_3 = np.zeros(3)
    
    for p in range(3):
        r = random()
        l_3[p] = -lbd * log(1 - r)
        r = random()
        fi_3[p] = 2 * pi * r
    
    # направляющие косинусы
    a_3 = np.zeros(3)
    b_3 = np.zeros(3)
    c_3 = np.zeros(3)
    
    for p in range(3):
        a_3[p] = sqrt(1 - square(cos_t_3[p])) * cos(fi_3[p])
        b_3[p] = sqrt(1 - square(cos_t_3[p])) * sin(fi_3[p])
        c_3[p] = cos_t_3[p]
    
    # координаты электрона
    x_3 = np.zeros(3)
    y_3 = np.zeros(3)
    z_3 = np.zeros(3)
    
    t_3 = np.zeros(3)
    
    for p in range(3):
        l_square_array[p][i][0] = square(x_3[p]) + square(y_3[p]) + square(z_3[p])
        t_array[p][0] = t_3[p]
        
        # заполнение координат в первый момент
        x_array[p][i][0] = 0
        y_array[p][i][0] = 0
        z_array[p][i][0] = 0
    
    # полёт i-ого электрона
    for j in range(1, n):
        a_ = np.zeros(3)
        b_ = np.zeros(3)
        c_ = np.zeros(3)
        
        for p in range(3):
            a_[p] = a_3[p]
            b_[p] = b_3[p]
            c_[p] = c_3[p]
        
            r = random()
            l_3[p] = -lbd * log(1 - r)
            
            r = random()
            fi_3[p] = 2 * pi * r
        
        r = random()
        cos_t_3[0] = 1 - 2 * r
        r = random()
        cos_t_3[1] = 1 - 2 * pow(r, 2/3)
        r = random()
        cos_t_3[2] = 2 * pow(1 - r, 2/3) - 1
        
        for p in range(3):
            if c_[p] != 1:
                sqrt_ = sqrt((1 - square(cos_t_3[p])) / (1 - square(c_[p])))
                
                a_3[p] = cos_t_3[p] * a_[p] - (b_[p] * sin(fi_3[p]) - 
                                               a_[p] * c_[p] * cos(fi_3[p])) * sqrt_
                b_3[p] = cos_t_3[p] * b_[p] + (a_[p] * sin(fi_3[p]) + 
                                               b_[p] * c_[p] * cos(fi_3[p])) * sqrt_
                c_3[p] = cos_t_3[p] * c_[p] - (1 - square(c_[p])) * cos(fi_3[p]) * sqrt_
                
            else:
                sin_t = sqrt(1 - square(cos_t[p]))
                
                a_3[p] = sin_t * cos(fi_3[p])
                b_3[p] = sin_t * sin(fi_3[p])
                c_3[p] = cos_t_3[p]
        
            x_3[p] += a_3[p] * l_3[p]
            y_3[p] += b_3[p] * l_3[p]
            z_3[p] += c_3[p] * l_3[p]
        
            t_3[p] += l_3[p] / v
        
            l_square_array[p][i][j] = square(x_3[p]) + square(y_3[p]) + square(z_3[p])
            t_array[p][j] = t_3[p]
            
            x_array[p][i][j] = x_3[p]
            y_array[p][i][j] = y_3[p]
            z_array[p][i][j] = z_3[p]

# нахождение усреднённого квадрата смещения K электронов для каждого столкновения
for p in range(3):
    average_l_square_3 = np.zeros(3)
    
    for j in range(n):
        
        for i in range(k):
            average_l_square_3[p] += l_square_array[p][i][j]
        
        average_l_square_3[p] /= k
        average_l_square_array[p][j] = average_l_square_3[p]

# массив среднеквадратичных смещений электронов
l_mean_square_array = np.zeros((3, n))

for p in range(3):
    for j in range(n):
        l_mean_square_array[p][j] = sqrt(average_l_square_array[p][j])

# теоретически рассчитанный квадрат смещения (по закону Эйнштейна)
l_theor_square_array = 2 * lbd * v * t_array[0]

# теоретически рассчитанное смещение (по закону Эйнштейна)
l_theor_array = sqrt(l_theor_square_array)

# построение графиков для K электронов
plt.figure(figsize = (12, 6))
plt.plot(t_array[0], l_mean_square_array[0], label = "σ(θ) = const")
plt.plot(t_array[1], l_mean_square_array[1], label = "σ(θ) = sin(" + r"$\frac{θ}{2}$" + ")")
plt.plot(t_array[2], l_mean_square_array[2], label = "σ(θ) = cos(" + r"$\frac{θ}{2}$" + ")")
plt.plot(t_array[0], l_theor_array, label = "теоретич. знач.")
plt.legend()
plt.title("Зависимость среднеквадратичного смещения от времени, K = " + str(k), fontsize = 18)
plt.xlabel("t, в единицах λ/v", fontsize = 14)
plt.ylabel("l, в единицах λ", fontsize = 14)
plt.show()

plt.figure(figsize = (12, 6))
plt.plot(t_array[0], average_l_square_array[0], label = "σ(θ) = const")
plt.plot(t_array[1], average_l_square_array[1], label = "σ(θ) = sin(" + r"$\frac{θ}{2}$" + ")")
plt.plot(t_array[2], average_l_square_array[2], label = "σ(θ) = cos(" + r"$\frac{θ}{2}$" + ")")
plt.plot(t_array[0], l_theor_square_array, label = "теоретич. знач.")
plt.legend()
plt.title("Зависимость среднего квадрата смещения от времени, K = " + str(k), fontsize = 18)
plt.xlabel("t, в единицах λ/v", fontsize = 14)
plt.ylabel("l², в единицах λ²", fontsize = 14)
plt.show()

last_moments = 3
step = 10

last_x_array = np.zeros((3, last_moments, k))
last_y_array = np.zeros((3, last_moments, k))
last_z_array = np.zeros((3, last_moments, k))

# Двумя строками ниже закомментирована команда, выводящая графики в отдельное окно
# Возможны ошибки в случае раскомментирования следующей строки

#%matplotlib qt

# графики для нескольких моментов времени
for p in range(3):
    for j in range(last_moments):
        for i in range(k):
            last_x_array[p][j][i] = x_array[p][i][n - 1 - step * (last_moments - j)]
            last_y_array[p][j][i] = y_array[p][i][n - 1 - step * (last_moments - j)]
            last_z_array[p][j][i] = z_array[p][i][n - 1 - step * (last_moments - j)]

num_p = 0

# for j in range(last_moments):
#     ax = Axes3D(plt.figure())
#     ax.scatter(last_x_array[num_p][j], last_y_array[num_p][j], last_z_array[num_p][j])
#     ax.scatter(0, 0, 0)
#     plt.title("Распределение электронов в пространстве, n = %d, k = %d" % (n, k))
#     plt.xlabel("x")
#     plt.ylabel("y")
#     plt.show()

# Распределение электронов по плоскости xOy

last_max_x = last_min_x = 0
last_max_y = last_min_y = 0

# Нахождение минимальных и максимальных координат (на xOy)
for i in range(k):
    if last_x_array[num_p][last_moments - 1][i] > last_max_x:
        last_max_x = last_x_array[num_p][last_moments - 1][i]
    elif last_x_array[num_p][last_moments - 1][i] < last_min_x:
        last_min_x = last_x_array[num_p][last_moments - 1][i]
    
    if last_y_array[num_p][last_moments - 1][i] > last_max_y:
        last_max_y = last_y_array[num_p][last_moments - 1][i]
    elif last_y_array[num_p][last_moments - 1][i] < last_min_y:
        last_min_y = last_y_array[num_p][last_moments - 1][i]

if last_max_x > last_max_y:
    last_max_y = last_max_x
else:
    last_max_x = last_max_y

if last_min_x < last_min_y:
    last_min_y = last_min_x
else:
    last_min_x = last_min_y

if last_max_x > - last_min_x:
    last_min_x = - last_max_x
else:
    last_max_x = - last_min_x

if last_max_y > - last_min_y:
    last_min_y = - last_max_y
else:
    last_max_y = - last_min_y

# Разбивка плоскости xOy на ячейки, ниже их количество
space_2D_x_size = space_2D_y_size = 20

# Физический размер ячейки (см. след. комментарий)
size_cell_x = (last_max_x - last_min_x) / (space_2D_x_size - 1)
size_cell_y = (last_max_y - last_min_y) / (space_2D_y_size - 1)

# Уточнение минимальных и максимальных координат (сдвиг на пол-ячейки от крайних электронов)
last_min_x -= size_cell_x / 2
last_max_x += size_cell_x / 2

last_min_y -= size_cell_y / 2
last_max_y += size_cell_y / 2

# Двумерная сетка, в каждой ячейке которой находится сколько-то электронов
space_2D = np.zeros((space_2D_x_size, space_2D_y_size))

# Определение номеров ячейки, в которой находится рассматриваемый электрон,
# добавление этого электрона в ячейку
for i in range(k):
    num_cell_x = int(np.floor((last_x_array[num_p][last_moments - 1][i] - last_min_x) / size_cell_x))
    num_cell_y = int(np.floor((last_y_array[num_p][last_moments - 1][i] - last_min_y) / size_cell_y))
    
    space_2D[num_cell_x][num_cell_y] += 1

# Двумерные массивы, необходимые для метода Axes3D.plot_surface()
x = np.linspace(last_min_x, last_max_x, space_2D_x_size)
y = np.linspace(last_min_y, last_max_y, space_2D_y_size)
x_2D, y_2D = np.meshgrid(x, y)

ax = Axes3D(plt.figure(figsize = (12, 7)))
ax.plot_surface(X = x_2D, Y = y_2D, Z = space_2D)

num_diff_x_in_cell = 10
dx = size_cell_x / num_diff_x_in_cell

num_diff_y_in_cell = 10
dy = size_cell_y / num_diff_y_in_cell

num_diff_z_in_volume = num_diff_x_in_cell * space_2D_x_size
dz = (last_max_x - last_min_x) / num_diff_z_in_volume

d = lbd * v / 3
t = t_array[num_p][n - 1]

space_2D_theor = np.zeros((space_2D_x_size, space_2D_y_size))

last_max_r = last_max_x
size_cell_r = size_cell_x
space_2D_r_size = int(space_2D_x_size / 2 + 1)
num_electrons_r = np.zeros(space_2D_r_size)

for idr in range(space_2D_r_size):
    for idx in range(num_diff_x_in_cell):
        x = idx * dx + idr * size_cell_r
        for idy in range(num_diff_y_in_cell):
            y = - size_cell_r / 2 + idy * dy
            for idz in range(num_diff_z_in_volume):
                z = - last_max_r + idz * dz
                _4dt = 4 * d * t
                num_electrons_r[idr] += dx * dy * k / ((pi * _4dt) ** 1.5) * exp(- (square(x) + 
                                                        square(y) + square(z)) / _4dt)

for idx in range(space_2D_x_size):
    for idy in range(space_2D_y_size):
        x = last_min_x + idx * size_cell_x
        y = last_min_y + idy * size_cell_y
        r = sqrt(square(x) + square(y))
        if r > last_max_r:
            space_2D_theor[idx][idy] = 0
        else:
            num_r = int(r / last_max_r * space_2D_r_size) - 1
            space_2D_theor[idx][idy] = num_electrons_r[num_r]
space_2D_theor[int(space_2D_x_size / 2)][int(space_2D_y_size / 2)] = num_electrons_r[0]

ax_2 = Axes3D(plt.figure(figsize = (12, 7)))
ax_2.plot_surface(X = x_2D, Y = y_2D, Z = space_2D_theor)

# for idx in range(space_2D_x_size):
#     for idy in range(space_2D_y_size):
#         for idx_2 in range(num_diff_x_in_cell):
#             x = last_min_x + idx * size_cell_x + idx_2 * dx
#             for idy_2 in range(num_diff_y_in_cell):
#                 y = last_min_y + idy * size_cell_y + idy_2 * dy
#                 for idz_2 in range(num_diff_z_in_volume):
#                     z = last_min_x + idz_2 * dz
#                     _4dt = 4 * d * t
#                     space_2D_theor[idx][idy] += dx * dy * k / ((pi * _4dt) ** 1.5) * exp(- (square(x) + square(y) + square(z)) / _4dt)

# ax_2 = Axes3D(plt.figure(figsize = (12, 7)))
# ax_2.plot_surface(X = x_2D, Y = y_2D, Z = space_2D_theor)
