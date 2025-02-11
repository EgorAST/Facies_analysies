import numpy as np
from typing import Union
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


def find_best_match(fragment: np.ndarray, full_curve: np.ndarray) -> Union[int, float]:
    fragment_length = len(fragment)
    best_match_index = -1
    min_area = -np.inf

    # Нормализуем фрагмент
    fragment_norm = normalize_data(fragment)

    # Вычисляем производную фрагмента
    fragment_derivative = np.gradient(fragment_norm)
    segment_derivative_ = None

    for i in range(len(full_curve) - fragment_length + 1):
        current_segment = full_curve[i:i + fragment_length]
        current_segment_norm = normalize_data(current_segment)

        # Вычисляем производную текущего сегмента
        segment_derivative = np.gradient(current_segment_norm)

        # Вычисляем корреляцию между производными
        correlation = np.corrcoef(fragment_derivative, segment_derivative)[0, 1]

        # Сохраняем лучшее совпадение
        if correlation > min_area:
            print(f"correlation {correlation}")
            segment_derivative_ = segment_derivative
            min_area = correlation
            best_match_index = i

    """plt.figure(figsize=(10, 6))

    freq = np.arange(len(fragment_derivative))
    # График амплитудного спектра первой последовательности
    plt.plot(freq, fragment_derivative, label="Производная Эталона", color="blue", alpha=1)

    # График амплитудного спектра второй последовательности
    plt.plot(freq, segment_derivative_, label="Производная Объекта", color="red", alpha=1)
    plt.title(f"корреляция = {min_area}")
    plt.legend()
    plt.grid(True)
    plt.show()"""

    return best_match_index, min_area