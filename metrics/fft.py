import numpy as np
from typing import Union
from scipy.fft import fft
from scipy.spatial.distance import cosine
import matplotlib.pyplot as plt
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


def fourier_similarity(sequence1: np.ndarray, sequence2: np.ndarray) -> float:
    """
    Оценивает сходство двух кривых на основе анализа Фурье.

    Параметры:
    ----------
    sequence1 : np.ndarray
        Первая последовательность (кривая).
    sequence2 : np.ndarray
        Вторая последовательность (кривая).

    Возвращает:
    ----------
    float
        Мера сходства двух кривых на основе Фурье-анализа (1 - идентичные, 0 -完全不同).
    """
    if len(sequence1) != len(sequence2):
        raise ValueError("Последовательности должны иметь одинаковую длину.")

    # Применяем быстрое преобразование Фурье (FFT) к обеим последовательностям
    fft_sequence1 = fft(sequence1)
    fft_sequence2 = fft(sequence2)

    # Вычисляем амплитудные спектры (модули комплексных чисел)
    amplitude_spectrum1 = np.abs(fft_sequence1)
    amplitude_spectrum2 = np.abs(fft_sequence2)

    param_1 = amplitude_spectrum1.copy()
    param_2 = amplitude_spectrum2.copy()

    # Нормируем амплитудные спектры
    amplitude_spectrum1 /= np.linalg.norm(amplitude_spectrum1)
    amplitude_spectrum2 /= np.linalg.norm(amplitude_spectrum2)

    correlation_matrix = np.corrcoef(amplitude_spectrum1, amplitude_spectrum2)
    correlation = correlation_matrix[0, 1]

    return correlation, param_1, param_2
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
    min_area = - np.inf
    fragment_norm = normalize_data(fragment)


    for i in range(len(full_curve) - fragment_length + 1):
        current_segment = full_curve[i:i + fragment_length]
        current_segment_norm = normalize_data(current_segment)
        distance, amplitude_spectrum1, amplitude_spectrum2 = fourier_similarity(current_segment_norm, fragment_norm)
        if distance > min_area:
            print(F"distance {distance}, fft")
            min_area = distance
            best_match_index = i



    return best_match_index, min_area