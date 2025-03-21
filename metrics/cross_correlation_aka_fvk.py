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


def cross_correlation(sequence1: np.ndarray, sequence2: np.ndarray) -> np.ndarray:
    """
    Вычисляет функцию взаимной корреляции (ФВК) между двумя последовательностями.

    Параметры:
    ----------
    sequence1 : np.ndarray
        Первая последовательность (кривая).
    sequence2 : np.ndarray
        Вторая последовательность (кривая).

    Возвращает:
    ----------
    np.ndarray
        Массив значений ФВК для всех возможных временных сдвигов.
    """
    # Вычисляем ФВК
    ccf = np.correlate(sequence1, sequence2, mode='full')
    return ccf


def find_best_match(fragment: np.ndarray, full_curve: np.ndarray) -> Union[int, float]:
    fragment_length = len(fragment)
    best_match_index = -1
    max_correlation = -np.inf  # Ищем максимум корреляции

    # Нормализуем фрагмент
    fragment_norm = normalize_data(fragment)

    for i in range(len(full_curve) - fragment_length + 1):
        # Выделяем текущий сегмент из полной кривой
        current_segment = full_curve[i:i + fragment_length]
        current_segment_norm = normalize_data(current_segment)

        # Вычисляем взаимную корреляцию
        correlation = np.correlate(fragment_norm, current_segment_norm, mode='valid')

        # Нормируем корреляцию
        norm_correlation = correlation / (np.sqrt(np.sum(fragment_norm ** 2) * np.sum(current_segment_norm ** 2)))
        correlation_value = norm_correlation[0]

        # Обновляем лучший результат, если текущая корреляция больше
        if correlation_value > max_correlation:
            max_correlation = correlation_value
            best_match_index = i

    return best_match_index, max_correlation