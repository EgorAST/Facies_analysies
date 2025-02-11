import numpy as np
import cv2
from scipy.spatial.distance import euclidean
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


def calculate_hu_moments(sequence: np.ndarray) -> np.ndarray:
    """
    Вычисляет Hu моменты для кривой.

    Параметры:
    ----------
    sequence : np.ndarray
        Последовательность (кривая), представленная как массив точек (x, y).

    Возвращает:
    ----------
    np.ndarray
        Массив из 7 Hu моментов.
    """
    if not isinstance(sequence, np.ndarray):
        sequence = np.array(sequence)

    curve = sequence.reshape(-1, 1, 2).astype(np.float32)

    # Вычисляем моменты
    moments = cv2.moments(curve)

    # Вычисляем Hu моменты
    hu_moments = cv2.HuMoments(moments).flatten()

    return hu_moments


def find_best_match(fragment: np.ndarray, full_curve: np.ndarray) -> Union[int, float]:
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
    fragment_norm = [(x, i) for i, x in enumerate(fragment_norm)]

    for i in range(len(full_curve) - fragment_length + 1):
        current_segment = full_curve[i:i + fragment_length]
        current_segment_norm = normalize_data(current_segment)

        current_segment_norm = [(x, i) for i, x in enumerate(current_segment_norm)]

        hu_moments1 = calculate_hu_moments(current_segment_norm)
        hu_moments2 = calculate_hu_moments(fragment_norm)

        if np.isnan(hu_moments1).any() or np.isinf(hu_moments1).any() or \
                np.isnan(hu_moments2).any() or np.isinf(hu_moments2).any():
            # Обработка случая с некорректными значениями
            distance = np.inf
        else:
            norm1 = np.linalg.norm(hu_moments1)
            norm2 = np.linalg.norm(hu_moments2)
            if norm1 == 0 or norm2 == 0:
                distance = np.inf
            else:
                hu_moments1_normalized = hu_moments1 / norm1
                hu_moments2_normalized = hu_moments2 / norm2
                distance = euclidean(hu_moments1_normalized, hu_moments2_normalized)

        if distance < min_area:

            min_area = distance
            best_match_index = i

    return best_match_index, min_area