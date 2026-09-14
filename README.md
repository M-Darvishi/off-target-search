# Approximate DNA Sequence Search

A **C++/Python** library for approximate DNA sequence searching using **Sliding Window** and **K-mer** based methods.

## Features

- Approximate DNA sequence matching with **Hamming distance**
- **Sliding Window** and **K-mer** search methods
- **Array** and **Hash-based** K-mer indexes
- Serial and **OpenMP** parallel implementations
- **Python** interface using **pybind11**
- Performance evaluation based on **time**, **candidates**, and **memory usage**
- Cache support for reusing built indexes

## Methods

- Sliding Window
- K-mer
  - Naive
  - Array
  - Hash

## Benchmarks

The project includes experiments for:

- Thread count
- K value
- Cache usage
- Memory consumption

## Data

Experiments use the **E. coli** reference genome.

---

# جستجوی تقریبی توالی DNA

یک کتابخانه **C++/Python** برای جستجوی تقریبی توالی‌های DNA با استفاده از روش‌های **Sliding Window** و **K-mer**.

## قابلیت‌ها

- جستجوی تقریبی بر اساس **Hamming Distance**
- روش‌های **Sliding Window** و **K-mer**
- Index مبتنی بر **Array** و **Hash**
- پیاده‌سازی سریال و موازی با **OpenMP**
- رابط **Python** با استفاده از **pybind11**
- ارزیابی **زمان اجرا**، **تعداد Candidateها** و **مصرف حافظه**
- پشتیبانی از **Cache** برای استفاده مجدد از Index

## روش‌ها

- Sliding Window
- K-mer
  - Naive
  - Array
  - Hash

## آزمایش‌ها

آزمایش‌های پروژه شامل بررسی موارد زیر است:

- تعداد Thread
- مقدار K
- تأثیر Cache
- مصرف حافظه

## داده

آزمایش‌ها با استفاده از ژنوم مرجع **E. coli** انجام شده‌اند.


---

![C++](https://img.shields.io/badge/C%2B%2B-17-00599C)
![OpenMP](https://img.shields.io/badge/OpenMP-Parallel_Computing-0066CC)
![Python](https://img.shields.io/badge/Python-3.x-3776AB)
![pybind11](https://img.shields.io/badge/pybind11-Python%2FC%2B%2B-FFD43B)
