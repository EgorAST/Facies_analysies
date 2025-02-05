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


def calculate_mse(sequence1: np.ndarray, sequence2: np.ndarray) -> float:
    """
    Вычисляет среднеквадратичное отклонение (MSE) между двумя последовательностями.

    Параметры:
    ----------
    sequence1 : np.ndarray
        Первая последовательность (кривая).
    sequence2 : np.ndarray
        Вторая последовательность (кривая).

    Возвращает:
    ----------
    float
        Среднеквадратичное отклонение между двумя последовательностями.
    """
    if len(sequence1) != len(sequence2):
        raise ValueError("Последовательности должны иметь одинаковую длину.")

    # Вычисляем квадрат разности между соответствующими точками
    squared_diff = np.square(sequence1 - sequence2)

    # Возвращаем среднее значение квадратов разностей
    return np.mean(squared_diff)


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


    for i in range(len(full_curve) - fragment_length + 1):
        current_segment = full_curve[i:i + fragment_length]
        current_segment_norm = normalize_data(current_segment)
        distance = calculate_mse(current_segment_norm, fragment_norm)
        if distance < min_area:
            min_area = distance
            best_match_index = i

    return best_match_index, min_area