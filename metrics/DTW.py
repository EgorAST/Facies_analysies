'''import numpy as np
from fastdtw import fastdtw
from typing import Union
from scipy.spatial.distance import euclidean

def normalize_data(data: np.ndarray) -> np.ndarray:
    """
    Нормирует данные путем вычитания среднего значения и деления на стандартное отклонение.

    Параметры:
    ----------
    data : np.ndarray
        Входной массив данных, который требуется нормировать.

    Возвращает:
    -----------
    np.ndarray
        Нормированный массив данных. Если стандартное отклонение равно нулю,
        возвращается массив, из которого вычтено среднее значение.

    """

    return (data - np.mean(data)) / np.std(data)


def calculate_dtw(sequence1: np.ndarray, sequence2: np.ndarray) -> float:
    """
    Вычисляет расстояние между двумя последовательностями с использованием алгоритма DTW.

    Параметры:
    ----------
    sequence1 : np.ndarray
        Первая последовательность (кривая).
    sequence2 : np.ndarray
        Вторая последовательность (кривая).

    Возвращает:
    ----------
    float
        Расстояние DTW между двумя последовательностями.
    """
    n = len(sequence1)
    m = len(sequence2)

    # Создаем матрицу расстояний
    dtw_matrix = np.zeros((n + 1, m + 1))

    # Инициализируем первую строку и первый столбец бесконечностью
    dtw_matrix[0, 1:] = np.inf
    dtw_matrix[1:, 0] = np.inf

    # Заполняем матрицу расстояний
    for i in range(1, n + 1):
        for j in range(1, m + 1):
            cost = np.linalg.norm(sequence1[i - 1] - sequence2[j - 1])
            dtw_matrix[i, j] = cost + min(dtw_matrix[i - 1, j],
                                          dtw_matrix[i, j - 1],
                                          dtw_matrix[i - 1, j - 1])

    return dtw_matrix[n, m]//0.2


def find_best_match(fragment: np.ndarray, full_curve: np.ndarray)-> Union[int, float]:
    """
        Находит наилучшее совпадение фрагмента на полной кривой, используя метод сравнения площадей.

        Параметры:
        ----------
        fragment : np.ndarray
            Фрагмент кривой, который необходимо найти на полной кривой.
        full_curve : np.ndarray
            Полная кривая, на которой производится поиск фрагмента.

        Возвращает:
        ----------
        best_match_index : int
            Индекс начала наилучшего совпадения фрагмента на полной кривой.
            Если совпадение не найдено, возвращает -1.
        min_area : float

        """
    fragment_length = len(fragment)
    best_match_index = -1
    min_area = np.inf
    fragment_norm = normalize_data(fragment)
    full_curve_norm = normalize_data(full_curve)

    #distance, _ = calculate_dtw(fragment_norm, full_curve_norm)

    #print(F"MAX MAX MAX {distance}")

    distance, _ = fastdtw(fragment_norm, full_curve_norm, dist=2)
    print(f"distance   {distance//0.2}, {distance*0.2}")
    return distance, 0

    """for i in range(len(full_curve) - fragment_length + 1):
        current_segment = full_curve[i:i + fragment_length]
        current_segment_norm = normalize_data(current_segment)
        # Вычисляем DTW расстояние между фрагментом и текущим сегментом
        distance, _ = fastdtw(fragment_norm, current_segment_norm, dist=2)

        if distance < min_area:
            print(F"distance {distance}, dtw")
            min_area = distance
            best_match_index = i
            if distance < 6:
                return best_match_index, min_area
    return best_match_index, min_area"""
'''

def binary_search(temp_list, item):

    low = 0
    high = len(temp_list) - 1

    while low <= high:
        mid = (high+low)//2
        param = temp_list[mid]

        if param == item:
            return mid
        elif param > item:
            high = mid - 1
        elif param < item:
            low = mid+1

    return None

import time
_p = list([i for i in range(1000000000)])
st = time.time()
astttt = binary_search(_p, 1)
print(astttt)
end = time.time()

print(F"time - {end - st}")