# Задаємо значення з вашої матриці помилок (Confusion Matrix)
tp = 54  # Імуногенні з високим Foreignness (>0.5)
fn = 29  # Імуногенні з низьким Foreignness (<=0.5)
fp = 33  # Неімуногенні з високим Foreignness (>0.5)
tn = 221  # Неімуногенні з низьким Foreignness (<=0.5)

total = tp + tn + fp + fn

# Розрахунок метрик
accuracy = (tp + tn) / total
precision = tp / (tp + fp)
recall = tp / (tp + fn)
f1_score = 2 * (precision * recall) / (precision + recall)

# Виведення результатів (у відсотках для зручності)
print(f"Загальна кількість пептидів: {total}")
print(f"Accuracy (Точність): {accuracy:.4f} ({accuracy * 100:.2f}%)")
print(f"Precision (Точність прогнозу): {precision:.4f} ({precision * 100:.2f}%)")
print(f"Recall (Повнота / Чутливість): {recall:.4f} ({recall * 100:.2f}%)")
print(f"F1 Score (Гармонійне середнє): {f1_score:.4f} ({f1_score * 100:.2f}%)")