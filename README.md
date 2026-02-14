
<div dir="rtl" style="direction: rtl; text-align: right;">

[![DOI](https://zenodo.org/badge/967615825.svg)](https://doi.org/10.5281/zenodo.18393840)

# M5: العلوم المفتوحة لظروف المحاصيل

> كيف يمكن استخدام مجموعات بيانات NASA لرسم خرائط ظروف المحاصيل؟

تركّز الوحدة الخامسة من [منهجنا المفتوح لعلوم المناخ](https://openclimatescience.github.io/curriculum) على كيفية بدء مشروع علمي حاسوبي قابل لإعادة الإنتاج، باستخدام ظروف المحاصيل كمثال موضوعي. **في نهاية هذه الوحدة، ينبغي أن تكون قادرًا على:**

- الوصول إلى مجموعات بيانات قائمة على الأقمار الصناعية حول إنتاجية النبات، حالته، والنتح‑التبخر، واستخدامها.
- حساب مؤشر رضا متطلبات المياه للمحاصيل الزراعية.


## المحتويات

1. [إنشاء خطة مشروع](https://github.com/OpenClimateScience/M5-Open-Science-for-Crops/blob/main/notebooks/01_Creating_a_Project_Plan.ipynb)
2. [إنشاء بيئة بحثية قابلة لإعادة الإنتاج](https://github.com/OpenClimateScience/M5-Open-Science-for-Crops/blob/main/notebooks/02_Creating_a_Reproducible_Research_Environment.ipynb)
3. [إعداد سير عمل معالجة البيانات](https://github.com/OpenClimateScience/M5-Open-Science-for-Crops/blob/main/notebooks/03_Setting_Up_Data_Processing_Workflows.ipynb)
4. [كتابة خوارزمية شفافة](https://github.com/OpenClimateScience/M5-Open-Science-for-Crops/blob/main/notebooks/04_Writing_a_Transparent_Algorithm.ipynb)


## البدء

[راجع دليل التثبيت هنا.](https://github.com/OpenClimateScience/M1-Open-Climate-Data/blob/main/HOW_TO_INSTALL.md)

يمكنك تشغيل دفاتر Jupyter في هذا المستودع باستخدام [Github Codespaces](https://docs.github.com/en/codespaces/overview) أو كحاوية تطوير في [VSCode Dev Container](https://code.visualstudio.com/docs/devcontainers/containers). بعد تشغيل الحاوية، شغّل Jupyter Notebook كما يلي:

```sh
# Create your own password when prompted
jupyter server password

# Then, launch Jupyter Notebook; enter your password when prompted
jupyter notebook
```

**يمكن تثبيت مكتبات Python المطلوبة للتمارين باستخدام مدير الحزم `pip`:**

```sh
pip install xarray netcdf4 dask
```


## نواتج التعلم

يغطي هذا المقرر الكفاءات الأساسية التالية في [علم البيانات الحاسوبي](https://github.com/OpenClimateScience/Core-Competencies/blob/main/ScienceCore-Competencies.md):

- توثيق العلاقة بين الشيفرة، والنتائج، والبيانات الوصفية (CC1.5)
- استخدام مدير حزم لتثبيت وإدارة الاعتماديات البرمجية (CC1.10)
- القدرة على توسيع نطاق سير عمل حاسوبي (CC2.6)
- فهم إصدارات البرمجيات وإدارة النسخ (CC4.4)

**بالإضافة إلى ذلك، سيتعلم المشاركون كيفية:**

- استخدام [Pixi](https://pixi.sh/latest/) لإدارة بيئة البرمجيات والاعتماديات.
- استخدام [Snakemake](https://snakemake.github.io/) لأتمتة المهام القابلة لإعادة الإنتاج.
- حساب مؤشر رضا متطلبات المياه (Water Requirements Satisfaction Index).


## شكر وتقدير

تم تمكين هذا المنهج بفضل منحة من برنامج NASA Transition to Open Science (TOPS) رقم 80NSSC23K0864، وهو جزء من [برنامج NASA للتحول إلى العلوم المفتوحة](https://nasa.github.io/Transform-to-Open-Science/).

</div>
