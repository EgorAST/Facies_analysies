# -*- coding: utf-8 -*-
#Программа для загрузки LAS файлов и выполнения фациального анализа в автоматическом режиме

import pandas as pd
import os
import lasio
import tkintermapview as tkmap
import tkinter as tk
from tkinter import ttk
from tkinter import filedialog
import matplotlib.pyplot as plt
from scipy.interpolate import griddata
import numpy as np
from metrics import (square_method, euclidean_distance, Manhattan_distance, Hausdorf_distance, DTW, RMS, MAE,
                     pearson_correlation, cross_correlation_aka_fvk, Koef_spearmanr, fft, hu_moments, derivative)
from functools import wraps
import time

MAX_DEPTH = 0
MIN_DEPTH = float("inf")
ZOOM_LVL = 1
NUM_COLUMN = 0

select_curve = False

curves = {}
well_value = {}
well_coord = {}


list_val = []
list_canvas = []
instances = []
list_canvas_scroll = []
scores = []
depths = []
list_obj_on_canvas_2 = []
rep_list = []
list_instance_val = []
abs_val = []
boundary_list = []

rep_name = ["G1t", "G1b", "G2 t", "G2 b", "k2", "К2 kamp", "k1al2", "k1al3+al2", "k1rep1", "k1rep2gl", "k1al1" ,
            "k1apt+vot", "k1apt+gl", "K1apt1-I", "K1apt1-II", "K1nk-III", "J-IV", "j3", "J2bj", "J2bjr1top", "J2bjr1bot", "J2baa"]

name_well = [1, 2, 3, 4, 5, 6,]


root = tk.Tk()
root.state('zoomed')
root.title("LAS reader")

menu = tk.Menu(root, font=("Inter", 16))
root.config(menu=menu)

notebook = ttk.Notebook(root)

tab1 = tk.Frame(notebook)
tab2 = tk.Frame(notebook)
tab3 = tk.Frame(notebook)

def determ_the_time_spent(func_0):
    @wraps(func_0)
    def func(*args, **kwargs):
        start_time = time.time()
        result = func_0(*args, **kwargs)
        end_time = time.time()
        print(F"Время затраченное на выполнение функции {func_0} = {end_time-start_time:.2f}")
        return result
    return func

class PS_curve():
    """
        Класс для визуализации и взаимодействия с кривыми данных в скважине.

        Атрибуты:
        ----------
        master : tkinter.Widget
            Родительский виджет, в котором будет отображаться кривая.
        depth : list
            Список, содержащий минимальную, максимальную глубину и шаг глубины.
        curve : dict
            Словарь, содержащий данные кривой, где ключ — название кривой, а значение — список значений.
        well_name : str
            Название скважины.
        set_curve : set
            Множество уникальных значений кривой.
        canvas_w : int
            Ширина канваса для отображения кривой.
        canvas_h : int
            Высота канваса для отображения кривой.
        scroll_value : int
            Текущее значение прокрутки канваса.
        list_point : list
            Список для хранения координат точек, выбранных на канвасе.
        list_curve : list
            Список значений кривой.
        delta_depth : float
            Шаг глубины, рассчитанный на основе высоты канваса.
        delta_val : float
            Шаг значений кривой, рассчитанный на основе ширины канваса.
        max_ : float
            Максимальное значение кривой.

        Методы:
        -------
        __init__(self, master, depth: list, curve: dict, well_name: str):
            Инициализирует объект PS_curve, настраивает интерфейс и отображает кривую.

        draw_grid(self):
            Рисует координатную сетку на канвасе.

        on_canvas_click(self, event):
            Обрабатывает событие клика на канвасе, сохраняя координаты точки.

        on_canvas_motion(self, event):
            Обрабатывает движение мыши с зажатой кнопкой, рисуя прямоугольник выделения.

        on_canvas_release(self, event):
            Обрабатывает отпускание кнопки мыши, вычисляет выделенный фрагмент кривой и сравнивает его с другими кривыми.

        move_main_canvas(self, event):
            Обрабатывает движение мыши по канвасу, обновляя отображение значений и глубины.

        scroll_canvas(self, event):
            Обрабатывает прокрутку колеса мыши, изменяя масштаб или прокручивая канвас.
        """

    global list_canvas_scroll


    def __init__(self, master, depth: list, curve: dict, well_name: str):
        """
            Инициализирует объект PS_curve.

            Параметры:
            ----------
            master : tkinter.Widget
                Родительский виджет, в котором будет отображаться кривая.
            depth : list
                Список, содержащий минимальную, максимальную глубину и шаг глубины.
            curve : dict
                Словарь, содержащий данные кривой, где ключ — название кривой, а значение — список значений.
            well_name : str
                Название скважины.
        """
        global MAX_DEPTH, MIN_DEPTH, canvas_depth, ZOOM_LVL, list_canvas, instances, canvas_2, list_canvas_scroll, scores

        self.scroll_value = 0
        instances.append(self)
        self.master = master
        self.depth = depth
        self.curve = curve
        self.well_name = well_name
        self.set_curve = set(list(self.curve.values())[0])

        self.canvas_w = 190
        self.canvas_h = self.master.winfo_height()


        if MAX_DEPTH < self.depth[1]:
            MAX_DEPTH = self.depth[1]
            Canvas_depth_.max_h = self.depth[1]

        if MIN_DEPTH > self.depth[0]:
            Canvas_depth_.min_h = self.depth[0]
            MIN_DEPTH = self.depth[0]

        self.canvas_ = tk.Canvas(self.master)
        self.canvas_.grid(column=NUM_COLUMN, row=0)

        self.name_frame = tk.Frame(self.canvas_, background="white", border=2)
        self.name_frame.pack(side=tk.TOP, fill=tk.BOTH, padx=1, pady=1)

        self.lbl_name_well = tk.Label(self.name_frame, text=f"   {self.well_name}   ")
        self.lbl_name_well.pack(fill=tk.BOTH)

        self.val_frame = tk.Frame(self.canvas_, background="white", border=2)
        self.val_frame.pack(fill=tk.BOTH, padx=1, pady=1)

        self.lbl_val_min = tk.Label(self.val_frame, text=f"{round(min(self.set_curve), 1)}", background='#F0F0F0')
        self.lbl_val_min.pack(fill=tk.BOTH, side=tk.LEFT)

        self.lbl_val_id = tk.Label(self.val_frame, text=f" - mV - ", background='#F0F0F0')
        self.lbl_val_id.pack(fill=tk.BOTH, expand=1, side=tk.LEFT)

        self.lbl_val_max = tk.Label(self.val_frame, text=f"{round(max(self.set_curve), 1)}", background='#F0F0F0')
        self.lbl_val_max.pack(fill=tk.BOTH, side=tk.RIGHT)

        self.main_frame = tk.Frame(self.canvas_, background="white", border=2)
        self.main_frame.pack(fill=tk.BOTH, padx=1, expand=1, pady=1)

        self.y_scrollbar = tk.Scrollbar(self.main_frame, orient='vertical')
        self.y_scrollbar.pack(side='right', fill='y')

        self.main_canvas = tk.Canvas(self.main_frame, bg="white", yscrollcommand=self.y_scrollbar.set, width=self.canvas_w,
                                     height=self.canvas_h)
        self.main_canvas.pack(side='left', expand=True)

        canvas_2.update_idletasks()
        canvas_2.config(scrollregion=canvas_2.bbox("all"))
        self.draw_grid()

        list_canvas.append(self.main_canvas)

        self.y_scrollbar.config(command=self.main_canvas.yview)

        self.main_canvas.bind("<Motion>", self.move_main_canvas)
        self.main_canvas.bind("<MouseWheel>", self.scroll_canvas)

        self.main_canvas.bind("<Button-1>", self.on_canvas_click)
        self.main_canvas.bind("<B1-Motion>", self.on_canvas_motion)
        self.main_canvas.bind("<ButtonRelease-1>", self.on_canvas_release)

        self.master.update()
        self.master.update_idletasks()

        self.main_canvas.update()
        self.main_canvas.update_idletasks()

        self.ps_point = list(curve.values())[0]

        self.delta_depth = (self.depth[1] - self.depth[0]) // self.canvas_h
        self.delta_val = self.canvas_w // (max(self.set_curve) - min(self.set_curve))

        self.list_point = []

        self.list_curve = list(curve.values())[0]

        len_list = len(self.ps_point) - 1
        self.max_ = 0

        for i, point in enumerate(self.ps_point):
            if i < len_list:
                x_1 = point
                y_1 = self.depth[0] + (.2 * i)
                x_2 = self.ps_point[i + 1]
                y_2 = self.depth[0] + (.2 * (i + 1))

                self.max_ = x_1 if self.max_ < x_1 else self.max_

                self.max_ = x_2 if self.max_ < x_2 else self.max_

                self.main_canvas.create_line(x_1, y_1, x_2, y_2, fill="red", tags="line", width=2)

    def draw_rep_on_curve(self):
        id = name_well.index(int(self.well_name))
        list_temp = list(rep_list[id].value())

        for y in list_temp:
            self.main_canvas.create_line(0, y, self.canvas_w, y,
                                         fill='blue',
                                         width=1, tags="coord_line")
    def draw_grid(self):
        """
        Наносит на канвас координатную сетку
        """

        for x in range(0, self.canvas_w, 15):
            self.main_canvas.create_line(x, 0, x, 3000,
                                         fill='gray',
                                         width=1, tags="coord_line")

        # Рисуем горизонтальные линии
        for y in range(0, 3000, 15):
            self.main_canvas.create_line(0, y, self.canvas_w, y,
                                         fill='gray',
                                         width=1, tags="coord_line")

    def on_canvas_click(self, event):
        if select_curve:
            x, y = self.main_canvas.canvasx(event.x), self.main_canvas.canvasy(event.y)
            self.list_point = [x, y]

    def on_canvas_motion(self, event):
        self.main_canvas.delete("rec")
        if select_curve:
            x, y = self.main_canvas.canvasx(event.x), self.main_canvas.canvasy(event.y)
            self.main_canvas.create_rectangle(self.list_point[0], self.list_point[1], x, y, outline='blue', fill='blue',
                                              stipple='gray50', tags="rec")
    @determ_the_time_spent
    def on_canvas_release(self, event):
        """
        Функция принимает top и bottom отметки глубины
        должна проходится по другим кривым и сранивать схожести фрагментов
        return: dict = {self.name_well:[facie_top, facie_bot]}
        """
        global num, scores

        x, y = self.main_canvas.canvasx(event.x), self.main_canvas.canvasy(event.y)
        self.main_canvas.delete("rec")
        point_bot, point_top = y, self.list_point[1]

        firs_val_slice = int((point_top - self.depth[0]) / self.depth[2])
        sec_val_slice = int((point_bot - self.depth[0]) / self.depth[2])
        self.list_point = []

        curve_data = self.ps_point[firs_val_slice:sec_val_slice]

        for instance in instances:
            score, plp = cross_correlation_aka_fvk.find_best_match(curve_data, instance.ps_point)
            len_frag = point_bot - point_top
            val = instance.depth[0] + (score * .2)
            instance.main_canvas.create_rectangle(0, val, 370, val+len_frag, outline='green', fill='green',
                                              stipple='gray50', tags="rec")


    def move_main_canvas(self, event):
        x, y = self.main_canvas.canvasx(event.x), self.main_canvas.canvasy(event.y)

        self.lbl_val_id.configure(text=f"{float(x)} mV", foreground="red")
        Depth_curve.draw_depth(Canvas_depth_, y_val = float(y), y_coord=int(event.y))

        if select_curve:
            self.main_canvas.configure(cursor="crosshair")
        else:
            self.main_canvas.configure(cursor="")

    def scroll_canvas(self, event):
        global ZOOM_LVL

        if event.state & 0x0004:
            if event.delta > 0:
                ZOOM_LVL *= 1.2
            else:
                ZOOM_LVL *= 0.8

            for instance in instances:
                instance.main_canvas.delete("line")
                for i, point in enumerate(instance.ps_point):
                    if i < len(instance.ps_point) - 1:
                        x_1 = point
                        y_1 = instance.depth[0] + (.2 * i) * ZOOM_LVL
                        x_2 = instance.ps_point[i + 1]
                        y_2 = instance.depth[0] + (.2 * (i + 1)) * ZOOM_LVL

                        instance.max_ = x_1 if instance.max_ < x_1 else instance.max_

                        instance.max_ = x_2 if instance.max_ < x_2 else instance.max_

                        instance.main_canvas.create_line(x_1, y_1, x_2, y_2, fill="red", tags="line")

            draw_rep()

        elif event.delta > 0:
            for canvas_ in list_canvas:
                self.scroll_value -= 1
                canvas_.yview_scroll(-1, "units")

        elif event.delta < 0:
            for canvas_ in list_canvas:
                self.scroll_value += 1
                canvas_.yview_scroll(1, "units")


class Depth_curve:
    """
        Класс для канваса с глубинами.
        ...
        Атрибуты
        --------
        max_h : int
            Максимальная глубина
        min_h : int
            Минимальная глубина
        step : int
            Шаг дискретизации
        Методы
        ------
        draw_depth(y_temp, y_fact):
            Отрисовывает значение глубины на холсте.
    """
    def __init__(self, master, max_h, min_h, step):
        self.max_h, self.min_h, self.step, self.master = max_h, min_h, step, master

        self.canvas = tk.Canvas(self.master, width=90)
        self.canvas.pack(side=tk.RIGHT, fill=tk.BOTH, padx=5, pady=5)

        self.name_frame = tk.Frame(self.canvas, background="white", border=2)
        self.name_frame.pack(side=tk.TOP, fill=tk.BOTH, padx=1, pady=1)

        self.lbl_name_well = tk.Label(self.name_frame, text=f"   DEPTH   ")
        self.lbl_name_well.pack(fill=tk.BOTH)

        self.val_frame = tk.Frame(self.canvas, background="white", border=2)
        self.val_frame.pack(fill=tk.BOTH, padx=1, pady=1)

        self.lbl_val_id = tk.Label(self.val_frame, text=f" m ", background='#F0F0F0')
        self.lbl_val_id.pack(fill=tk.BOTH, expand=1, side=tk.LEFT)

        self.master.update()
        self.master.update_idletasks()

    def draw_depth(self, y_val: float, y_coord: int):
        """
            Отрисовывает значение глубины на холсте.

            Удаляет предыдущие элементы с тегом "depth_val" и создает новые:
            - Текстовое значение глубины (y_temp) в указанной позиции.
            - Две линии по бокам от текстового значения.
        """

        self.canvas.delete("depth_val")

        self.canvas.create_text(33, y_coord + 52, text=str(y_val), tags="depth_val")
        self.canvas.create_line(4, y_coord + 52, 14, y_coord + 52, tags="depth_val")
        self.canvas.create_line(54, y_coord + 52, 64, y_coord + 52, tags="depth_val")


def open_las_file(full_path=None):
    """
        Открывает и загружает данные из LAS-файла.
        Если путь к файлу не указан, открывает диалоговое окно для выбора файла.
        Загружает данные из LAS-файла, извлекает кривые и глубины, сохраняет их в глобальные переменные
        и отображает кривые на холсте, создавая экземпляры класса PS_curve.

        Параметры:
        ----------
        full_path : str, optional
            Полный путь к LAS-файлу. Если не указан, открывается диалоговое окно для выбора файла.


        Глобальные переменные:
        ----------------------
        curves : dict
            Словарь, содержащий данные кривых (кроме DEPT).
        depth_canvas : object
            Объект холста для отображения глубин.
        depths : list
            Список, содержащий минимальную, максимальную глубину и шаг.
        list_obj_on_canvas_2 : list
            Список объектов на холсте.
        NUM_COLUMN : int
            Счетчик колонок.

        """
    global curves, depth_canvas, depths, list_obj_on_canvas_2, NUM_COLUMN

    if 1:
        if full_path is None:
            file_path = tk.filedialog.askopenfilename(initialdir=".",
                                                      filetypes=[("LAS File", ".las")])
        else:
            file_path = full_path

        if file_path:

            name_ = (file_path.split('\\')[-1]).split(".")[0]

            las = lasio.read(file_path)
            keys = las.keys()

            # Создание объектов CurveCanvas для каждой кривой, кроме DEPT
            for i, key in enumerate(keys):
                if key != 'DEPT':
                    curve_data = np.array(las[key].data)  # Преобразуем в NumPy массив
                    curves[key] = curve_data  # Сохраняем данные в глобальном словаре
                else:
                    dept_data = np.array(las[key].data)
                    min_d, max_d, step = min(dept_data), max(dept_data), round((dept_data[1] - dept_data[0]), 2)
                    depths = [min_d, max_d, step]

            PS_curve(canvas_2, depths, curves, name_)
            NUM_COLUMN += 1


def open_die_las_file():
    """
        Поочередная загрузка LAS файлов из выбранной директории.

        Функция открывает диалоговое окно для выбора директории, после чего загружает все файлы с расширением `.las`
        из этой директории. Каждый LAS файл открывается с помощью функции `open_las_file`.
    """
    try:
        file_path = tk.filedialog.askdirectory(initialdir=".")
        files = os.listdir(file_path)

        if files:
            for file in files:
                if file.lower().endswith(".las"):
                    full_path = os.path.join(file_path, file)
                    open_las_file(full_path)

    except:
        tk.messagebox.showinfo("Ошибка", f"При загрузке файла произошла неизвестная ошибка")


def load_coord_well():
    """
        Загружает координаты скважин из текстового файла или Excel-файла и отображает их на карте.

        Функция открывает диалоговое окно для выбора файла с координатами скважин. Поддерживаются файлы формата `.txt` и `.xlsx`.
        После выбора файла, данные из него считываются, и координаты скважин добавляются на карту с помощью `gmap_widget.set_marker`.
        Координаты также сохраняются в глобальной переменной `well_coord` в формате словаря, где ключ — название скважины, а значение — кортеж с широтой и долготой.
    """
    global well_coord
    try:
        file_path = tk.filedialog.askopenfilename(initialdir=".",
                                                  filetypes=[("text File", ".txt"), ("Excel File", ".xlsx")])
        if file_path:
            with open(file_path, "r") as file:
                for string in file:
                    if len(string) > 1:
                        string = string.split()
                        gmap_widget.set_marker(float(string[-2]), float(string[-1]), text=f"{string[0]}")
                        well_coord[string[0]] = (string[-2], string[-1])
    except:
        tk.messagebox.showinfo("Ошибка", f"При загрузке файла произошла неизвестная ошибка")
        return

    else:
        full_well_value(button_frame_3)



def full_well_value(frame):
    global well_value, list_instance_val, abs_val

    for instance in instances:
        list_instance_val.append(instance.depth[0])


    for i, well in enumerate(list(well_coord)):
        well_value[str(well)] = list_val[i] * .2 + list_instance_val[i]
        abs_val.append(list_val[i] * .2 + list_instance_val[i])

    x = []
    y = []
    z = []

    for key in well_coord.keys():
        lat, lon = well_coord[key]
        x.append(float(lon))
        y.append(float(lat))
        z.append(well_value[key])

    x = np.array(x)
    y = np.array(y)
    z = np.array(z)

    # Создаем сетку для интерполяции
    grid_x, grid_y = np.mgrid[min(x):max(x):100j, min(y):max(y):100j]

    # Интерполяция (линейная, кубическая или nearest)
    grid_z = griddata((x, y), z, (grid_x, grid_y), method='linear')

    # Построение карты
    plt.figure(figsize=(10, 6))
    contour = plt.contourf(grid_x, grid_y,  grid_z, levels=20,
                           cmap='terrain_r')  # Меняем grid_x и grid_y местами, если нужно

    plt.colorbar(contour, label='Глубина м.')
    plt.scatter(x, y, c='red', s=50, edgecolor='k', label='Скважины')  # Отображение исходных точек
    plt.xlabel('Широта (Y)')
    plt.ylabel('Долгота (X)')
    plt.title('Поверхность')
    plt.legend()
    plt.grid(True)
    plt.show()


def load_rep():
    global rep_list

    try:
        file_path = tk.filedialog.askopenfilename(initialdir=".",
                                                  filetypes=[("Excel File", ".xlsx")])

        if file_path:
            if file_path.endswith(".xlsx"):
                # Чтение файла Excel
                df = pd.read_excel(file_path)
                for index, row in df.iterrows():
                    if index > 0:
                        rep_dict = {}

                        for i, r in enumerate(row):
                            if isinstance(r, float):
                                rep_dict[rep_name[i-1]] = r
                        rep_list.append(rep_dict)

    except:
        tk.messagebox.showinfo("Ошибка", f"При загрузке файла произошла неизвестная ошибка")

    finally:
        draw_rep()

def draw_rep():
    for instace in instances:

        instace.main_canvas.delete("rep")
        name = int(instace.well_name.split("_")[0])
        id_curve = int(name_well.index(name))

        list_temp = list(rep_list[id_curve].values())
        for i, y in enumerate(list_temp):
            instace.main_canvas.create_line(0, y, instace.canvas_w, y,
                                         fill='green',
                                         width=1, tags="rep")

            instace.main_canvas.create_text(33, y, text=str(rep_name[i]),  tags="rep")


def load_bondary():
    """
        Загружает граничные точки из файла и устанавливает их на карте.

        Функция открывает диалоговое окно для выбора файла (текстового или Excel),
        считывает данные из файла и добавляет координаты граничных точек в глобальный
        список `boundary_list`. После этого функция устанавливает полигон на карте
        с использованием этих точек.

        Глобальные переменные:
            boundary_list (list): Список для хранения координат граничных точек.
    """
    global boundary_list

    try:
        file_path = tk.filedialog.askopenfilename(initialdir=".",
                                                  filetypes=[("text File", ".txt"), ("Excel File", ".xlsx")])
        if file_path:

            with open(file_path, "r") as file:

                for string in file:

                    if len(string) > 1:
                        string = string.split()
                        boundary_list.append([float(string[-2]), float(string[-1])])

                gmap_widget.set_polygon(boundary_list)

    except:
        tk.messagebox.showinfo("Ошибка", f"При загрузке файла произошла неизвестная ошибка")

def on_key_press(event):
    global select_curve
    if event.keysym == "1":
        select_curve = not select_curve


# Создаем меню "Файл"
file_menu = tk.Menu(menu, tearoff=0)
menu.add_cascade(label="Файл", menu=file_menu, font=("Inter", 10))
file_menu.add_command(label="Открыть LAS файл", command=open_las_file, font=("Inter", 10))
file_menu.add_command(label="Открыть директорию с LAS файлами", command=open_die_las_file, font=("Inter", 10))
file_menu.add_command(label="Загрузить координаты для скважин", command=load_coord_well, font=("Inter", 10))
file_menu.add_command(label="Загрузить граница лицензионного участка", command=load_bondary, font=("Inter", 10))
file_menu.add_command(label="Загрузить отбивки", command=load_rep, font=("Inter", 10))


notebook.add(tab1, text='Карта')

button_frame = tk.Frame(tab1)
button_frame.pack(side=tk.TOP, fill=tk.X)

gmap_widget = tkmap.TkinterMapView(tab1, corner_radius=15)
gmap_widget.pack(fill=tk.BOTH, padx=5, pady=(4, 1), expand=True)

# ___

notebook.add(tab2, text='Каротаж')

button_frame_2 = tk.Frame(tab2)
button_frame_2.pack(side=tk.TOP, fill=tk.X)

left_frame = tk.Frame(tab2)
left_frame.pack(fill=tk.BOTH, side=tk.LEFT)

right_frame = tk.Frame(tab2)
right_frame.pack(fill=tk.BOTH, side=tk.RIGHT, expand=1)

x_scrollbar_2 = tk.Scrollbar(tab2, orient='horizontal')
x_scrollbar_2.pack(side='bottom', fill='x')

y_scrollbar_2 = tk.Scrollbar(tab2, orient='vertical')
y_scrollbar_2.pack(side='right', fill='y', pady=(10, 0))

canvas_depth = tk.Canvas(left_frame, bg="#FFFFFF", width=100)
canvas_depth.pack(side=tk.LEFT, fill=tk.BOTH, padx=5, pady=(5, 21))

Canvas_depth_ = Depth_curve(canvas_depth, 1, 0, .2)

canvas_2 = tk.Canvas(tab2, bg="#FFFFFF", width=10000, height=1000)
canvas_2.pack(padx=5, pady=5)

canvas_2.config(scrollregion=canvas_2.bbox("all"))
canvas_2.bind("<KeyPress>", on_key_press)
notebook.add(tab3, text='Поверхность')

button_frame_3 = tk.Frame(tab3)
button_frame_3.pack(expand=1, fill=tk.BOTH)

notebook.pack(expand=1, fill='both')




if __name__ == "__main__":
    root.mainloop()