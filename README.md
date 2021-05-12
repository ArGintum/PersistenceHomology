# PersistenceHomology

## Зависимости:
Необходимо наличие пакетов Boost и TBB (http://www.github.com/oneapi-src/oneTBB )

## Сборка:
Запустить make из папки с файлами Makefile и PersistenceCalculator.cpp

## Запуск
```sh
$ ./PersistenceCalculator [Имя файла с входными данными] [Параметры]
```

Доступные Параметры:
-  "--output": Имя выходного файла (по умолчанию output.txt )
-  "--format": Формат входного файла:  matrix (матрица расстояний, по умолчанию) или cloud (облако точек в евклидовом пространстве)
-  "--dim <k>": Максимальная группа гомологий, которая будет рассчитана (по умолчанию 2)
- "--threshold <t>": Максимальный порог расстояния (по умолчанию рассичтывается из входных данных)
- "--prep": Принудительное использование EdgeCollapse (по умолчанию, необходимость этого определеяется программой)
- "--noprep": Принудительный запрет на использование EdgeCollapse (по умолчанию, необходимость этого определеяется программой)
- "--threads <n>": Количество потоков для распараллеливания (по умолчанию определяется программой)

