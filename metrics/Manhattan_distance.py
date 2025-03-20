import numpy as np
from typing import Union


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

    mean = np.mean(data)
    std = np.std(data)
    if std == 0:
        return data - mean
    return (data - mean) / std


def calculate_area(segment1: np.ndarray, segment2: np.ndarray) -> float | int:
    """
    Функция определяет сумму Манхэттенсих расстояний между точкам
    :param segment1: Фрагмент кривой
    :param segment2: Кривая
    :return:
    """
    return np.sum(abs(segment1 - segment2))

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
            Минимальная площадь между нормализованным фрагментом и соответствующим сегментом полной кривой.
            Чем меньше площадь, тем лучше совпадение.

        """
    fragment_length = len(fragment)
    best_match_index = -1
    min_area = np.inf
    fragment_norm = normalize_data(fragment)


    for i in range(len(full_curve) - fragment_length + 1):
        current_segment = full_curve[i:i + fragment_length]
        current_segment_norm = normalize_data(current_segment)
        distance = calculate_area(current_segment_norm, fragment_norm)
        #print(F"dist {distance}")
        if distance < min_area:
            min_area = distance
            best_match_index = i

    return best_match_index, min_area