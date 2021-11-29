# BLAST representation by Sergei Nabatov

Цель проекта - представление BLAST для FPGA с использованием High Level Synthesis by Vitis HLS

## Description

Данный проект представляет собой реиспользование алгоритма BLAST в следующей последовательности:

1) Индексация БД (в том числе и сортировка)
2) Поиск short alignments с заданной строкой с использованием бинарного дерева
3) Для оценки используем Smith waterman алгоритм

Со стороны софта подаются две строки с определенными длинами (на данный момент захардкожены): строка и БД. 

Результатом выполнения является матрица схожести между двумя строками.

## Getting started

### Dependencies

* Для симуляции использовалась Vitis HLS 2020.2.
* Компилятор gcc 9.3.0

### Existing problems

- [ ] Проблема с сортировкой: она не экспортируется в RTL (по крайней мере на моем ПК, ему не хватает мощности). **Возможное решение**: создать два RTL блока в Vivado, в первом будет проводиться индексация БД, а во втором уже непосредственно расчет BLAST. 
- [ ] Необходима оптимизация для SW алгоритм и для индексации относительно hw разработки
- [ ] Необходимо рассчитать максимальное кол-во элементов, которые можно будет принимать на стороне hw в массивы и структуры
- [ ] Имплементация GSP совместно с Smith Waterman?

## Authors

Contributors names and contact info

ex. Sergei Nabatov [serynabatov@gmail.com]

## Version History

* 0.1 
    * Intial Release

## License 

This project is licensed under the Apache License - see the LICENSE.md file for details

## Acknowledgments

Code snippets used:
* [necst] (https://github.com/necst/coursera-sdaccel-practice/blob/master/Week-4/1-The_Smith-Waterman_example_in_details/1-On_how_to_optimize_the_Smith-Waterman_solution/1-A_first_implementation/final_code/maincl.cpp)

The algorithm has been built using these lectures:
* [Bionformatics Resource Core by Universidad de Puerto Rico] https://www.youtube.com/watch?v=jzSIC2UzxZ4&ab_channel=PR-INBREBiRC%5BBioinformaticsResourcesCore%5D
* [Rob Edwards San Diego State University Computational Genomics Class] https://www.youtube.com/watch?v=2V9HNxbWUMg&list=PLpPXw4zFa0uLMHwSZ7DMeLGjIUgo1IBbn&index=171&ab_channel=RobEdwards
