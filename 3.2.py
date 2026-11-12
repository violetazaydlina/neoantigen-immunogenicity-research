hla_data = {
    "HLA-A*01:01": {"tp": 5, "fn": 1, "fp": 8, "tn": 77},
    "HLA-A*02:01": {"tp": 30, "fn": 11, "fp": 15, "tn": 65},
    "HLA-A*02:05": {"tp": 1, "fn": 0, "fp": 0, "tn": 0},
    "HLA-A*03:01": {"tp": 5, "fn": 4, "fp": 0, "tn": 0},
    "HLA-A*11:01": {"tp": 8, "fn": 3, "fp": 1, "tn": 40},
    "HLA-A*24:02": {"tp": 1, "fn": 2, "fp": 0, "tn": 1},
    "HLA-A*26:01": {"tp": 0, "fn": 0, "fp": 9, "tn": 38},
    "HLA-A*29:02": {"tp": 3, "fn": 3, "fp": 0, "tn": 0},
    "HLA-A*30:01": {"tp": 0, "fn": 1, "fp": 0, "tn": 0},
    "HLA-A*31:01": {"tp": 0, "fn": 2, "fp": 0, "tn": 0},
    "HLA-A*32:01": {"tp": 0, "fn": 1, "fp": 0, "tn": 0},
    "HLA-A*68:01": {"tp": 1, "fn": 0, "fp": 0, "tn": 0},
    "HLA-A*68:02": {"tp": 0, "fn": 1, "fp": 0, "tn": 0},
}

for hla, vals in hla_data.items():
    tp = vals["tp"]
    fn = vals["fn"]
    fp = vals["fp"]
    tn = vals["tn"]

    total = tp + tn + fp + fn
    if total == 0:
        continue

    # Розрахунок метрик (з перевіркою на нуль у знаменниках)
    accuracy = (tp + tn) / total
    precision = tp / (tp + fp) if (tp + fp) > 0 else 0
    recall = tp / (tp + fn) if (tp + fn) > 0 else 0
    f1_score = (
        (2 * precision * recall) / (precision + recall)
        if (precision + recall) > 0
        else 0
    )

    # Виведення результатів для кожної алелі
    print(f"\n[HLA: {hla}]")
    print(f"Загальна кількість пептидів: {total}")
    print(f"Accuracy (Точність): {accuracy:.4f} ({accuracy * 100:.2f}%)")
    print(f"Precision (Точність прогнозу): {precision:.4f} ({precision * 100:.2f}%)")
    print(f"Recall (Повнота / Чутливість): {recall:.4f} ({recall * 100:.2f}%)")
    print(f"F1 Score (Гармонійне середнє): {f1_score:.4f} ({f1_score * 100:.2f}%)")