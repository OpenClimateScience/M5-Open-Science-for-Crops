
<div dir="rtl" style="direction: rtl; text-align: right;">

<p dir="ltr" style="text-align: left;">
<a href="https://doi.org/10.5281/zenodo.18393840">
<img src="https://zenodo.org/badge/967615825.svg" alt="DOI">
</a>
</p>

<h1 style="direction: rtl; text-align: right;">
M5: العلوم المفتوحة لظروف المحاصيل
</h1>

<blockquote style="direction: rtl; text-align: right;">
كيف يمكن استخدام مجموعات بيانات NASA لرسم خرائط ظروف المحاصيل؟
</blockquote>

<p style="direction: rtl; text-align: right;">
تركّز الوحدة الخامسة من 
<a href="https://openclimatescience.github.io/curriculum">منهجنا المفتوح لعلوم المناخ</a>
على كيفية بدء مشروع علمي حاسوبي قابل لإعادة الإنتاج، باستخدام ظروف المحاصيل كمثال موضوعي.
<strong>في نهاية هذه الوحدة، ينبغي أن تكون قادرًا على:</strong>
</p>

<ul style="direction: rtl; text-align: right;">
<li>الوصول إلى مجموعات بيانات قائمة على الأقمار الصناعية حول إنتاجية النبات، حالته، والنتح-التبخر، واستخدامها.</li>
<li>حساب مؤشر رضا متطلبات المياه للمحاصيل الزراعية.</li>
</ul>

<h2 style="direction: rtl; text-align: right;">المحتويات</h2>

<ol style="direction: rtl; text-align: right;">
<li><a href="https://github.com/OpenClimateScience/M5-Open-Science-for-Crops/blob/main/notebooks/01_Creating_a_Project_Plan.ipynb">إنشاء خطة مشروع</a></li>
<li><a href="https://github.com/OpenClimateScience/M5-Open-Science-for-Crops/blob/main/notebooks/02_Creating_a_Reproducible_Research_Environment.ipynb">إنشاء بيئة بحثية قابلة لإعادة الإنتاج</a></li>
<li><a href="https://github.com/OpenClimateScience/M5-Open-Science-for-Crops/blob/main/notebooks/03_Setting_Up_Data_Processing_Workflows.ipynb">إعداد سير عمل معالجة البيانات</a></li>
<li><a href="https://github.com/OpenClimateScience/M5-Open-Science-for-Crops/blob/main/notebooks/04_Writing_a_Transparent_Algorithm.ipynb">كتابة خوارزمية شفافة</a></li>
</ol>

<h2 style="direction: rtl; text-align: right;">البدء</h2>

<p style="direction: rtl; text-align: right;">
<a href="https://github.com/OpenClimateScience/M1-Open-Climate-Data/blob/main/HOW_TO_INSTALL.md">
راجع دليل التثبيت هنا
</a>
</p>

<p style="direction: rtl; text-align: right;">
يمكنك تشغيل دفاتر Jupyter باستخدام 
<a href="https://docs.github.com/en/codespaces/overview">Github Codespaces</a>
أو كحاوية تطوير في 
<a href="https://code.visualstudio.com/docs/devcontainers/containers">VSCode Dev Container</a>.
</p>

<pre dir="ltr">
jupyter server password
jupyter notebook
</pre>

<p style="direction: rtl; text-align: right;">
يمكن تثبيت مكتبات Python المطلوبة باستخدام مدير الحزم <code>pip</code>:
</p>

<pre dir="ltr">
pip install xarray netcdf4 dask
</pre>

<h2 style="direction: rtl; text-align: right;">نواتج التعلم</h2>

<ul style="direction: rtl; text-align: right;">
<li>توثيق العلاقة بين الشيفرة، والنتائج، والبيانات الوصفية (CC1.5)</li>
<li>استخدام مدير حزم لإدارة الاعتماديات البرمجية (CC1.10)</li>
<li>القدرة على توسيع نطاق سير عمل حاسوبي (CC2.6)</li>
<li>فهم إصدارات البرمجيات وإدارة النسخ (CC4.4)</li>
</ul>

<p style="direction: rtl; text-align: right;">
بالإضافة إلى ذلك، سيتعلم المشاركون كيفية:
</p>

<ul style="direction: rtl; text-align: right;">
<li>استخدام <a href="https://pixi.sh/latest/">Pixi</a> لإدارة بيئة البرمجيات.</li>
<li>استخدام <a href="https://snakemake.github.io/">Snakemake</a> لأتمتة المهام.</li>
<li>حساب مؤشر رضا متطلبات المياه (Water Requirements Satisfaction Index).</li>
</ul>

<h2 style="direction: rtl; text-align: right;">شكر وتقدير</h2>

<p style="direction: rtl; text-align: right;">
تم تمكين هذا المنهج بفضل منحة من برنامج NASA Transition to Open Science (TOPS)
رقم 80NSSC23K0864، وهو جزء من 
<a href="https://nasa.github.io/Transform-to-Open-Science/">برنامج NASA للتحول إلى العلوم المفتوحة</a>.
</p>

</div>
