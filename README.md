## Графики зависимости выполнения умножения CSR и плотной матрицы на вектор от времени:
| График 1 | График 2 |
|----------|----------|
| ![graph](https://github.com/user-attachments/assets/66df29ae-c679-462f-8489-65519edeb6c0) | ![graph_log](https://github.com/user-attachments/assets/bf78f32c-b642-42e9-8be6-6c47459acb55)
Видно, что CSR матрица умножается на вектор быстрее, чем плотная, в случая когда она разрежена.

## График зависимости невязки от номера итерации при решении СЛАУ для матрицы 5 на 5
![iterate_method](https://github.com/user-attachments/assets/378ccd78-ff85-455b-9a13-f7cb63323e9c)

## График зависимости количества итераций для сходимости метода SOR от параметра w
![SOR_of_w](https://github.com/user-attachments/assets/2ddd070e-ab39-4eff-9b12-15ab2a0b0d3c)
Видно, что у грфика есть минимум, что как раз и соответствует оптимальному значению параметра w.

## График сравнения скорости сходимости по итерациям и по времени ускоренных и обычных методов для матрицы 5 на 5
![iterate_method_vs_accel](https://github.com/user-attachments/assets/0d8967d1-e728-4f2f-b086-5136d587b6e7)
Видно что ускоренные методы сходятся быстрее по итерациям, но не по времени.

## График зависимости невязки от номера итерации при решении СЛАУ для матрицы 2500 на 2500
![iterate_method_new](https://github.com/user-attachments/assets/a9a9a2c3-84c2-4d96-8c8e-98c787faed4b)

## Сравнение скорости сходимости BiCG и CGS
![BiCG_vs_CGS](https://github.com/user-attachments/assets/8991a783-ee9b-413b-8e6e-cf48abfcf14c)
