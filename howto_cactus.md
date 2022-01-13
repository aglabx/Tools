# Как поставить кактус

Идем в гугл набираем "progressive cactus". 

Находим https://github.com/glennhickey/progressiveCactus, переходим видим, что там написано:

```text
IMPORTANT: Progressive Cactus has moved here:
https://github.com/ComparativeGenomicsToolkit/cactus
This version a) is actively maintained and developed and b) supports cloud computing platforms by using Toil in place of JobTree
```

Переходим.

Идем в секцию установка. Там видим

```text
Installation Overview
There are many different ways to install and run Cactus:

Docker Image
Precompiled Binaries
Build From Source
Python Install with Docker Binaries
```

Докер нам не подходит по религиозным соображениям. Для сборки из исходников нам будет нужен судо, его у нас нет. Поэтому наш путь бинарники.

Переходим по ссылке https://github.com/ComparativeGenomicsToolkit/cactus#precompiled-binaries

Видим:

```text
Precompiled Binaries
Precompiled binaries can be found on the Releases Page. Download by clicking the cactus-bin-vx.x.x.tar.gz and install following the instructions in the BIN-INSTALL.md link.
```

Жамкаем Releases Page. Это ссылка https://github.com/ComparativeGenomicsToolkit/cactus/releases

У нас линукс и выбираем Pre-compiled Binaries Linux Tarball

Наша ссылка для скачивания https://github.com/ComparativeGenomicsToolkit/cactus/releases/download/v2.0.4/cactus-bin-v2.0.4.tar.gz

Наша инструкция тут же Install instructions in BIN-INSTALL.md https://github.com/ComparativeGenomicsToolkit/cactus/blob/v2.0.4/BIN-INSTALL.md

Переходим.

Идем на сервер.

**Выходим из конды и других окружений**

**ИЗ ВСЕХ. ОНИ МОГУТ БЫТЬ ВЛОЖЕННЫМИ**

Слева в баше не должно быть имен окружений.

Дальше ставим ровно как написано в инструкции

**Но должно быть virtualenv. Если его нет то будет ошибка тогда попросить администратора поставить sudo apt install virtualenv**

**Предварительно лучше проверить версию rdflib-jsonld. Важно чтобы rdflib-jsonld<0.6.0 rdflib-jsonld>=0.3.0. Если версия именно от 0.6.0, то необходимо ГЛОБАЛЬНО установить rdflib-jsonld==0.5.0**

Устанавливаем, если необходимо, rdflib-jsonld:
```bash
pip install rdflib-jsonld==0.5.0
```
Дальше уже собираем кактус

```bash
wget https://github.com/ComparativeGenomicsToolkit/cactus/releases/download/v2.0.4/cactus-bin-v2.0.4.tar.gz
tar -xzf "cactus-bin-v2.0.4.tar.gz"
cd "cactus-bin-v2.0.4.tar.gz"

virtualenv -p python3.6 venv
source venv/bin/activate
pip install -U setuptools pip
pip install -U -r ./toil-requirement.txt
pip install -U .
export PATH=$(pwd)/bin:$PATH
export PYTHONPATH=$(pwd)/lib:$PYTHONPATH

cd bin && for i in wigToBigWig faToTwoBit bedToBigBed bigBedToBed bedSort hgGcPercent; do wget -q http://hgdownload.cse.ucsc.edu/admin/exe/linux.x86_64/${i}; chmod ugo+x ${i}; done
```

В процессе оно чуть будет сыпать ошибки. Ща поправим.

Дальше ставим сам кактус.

```bash
python setup.py install
pip install BioPython
pip install cigar
pip install cython
pip install pytest
```

Тестим работает ли

```bash
cactus ./jobstore ./examples/evolverMammals.txt ./evolverMammals.hal --realTimeLogging
```

Конец.

Для входа в окружения надо:

```bash
cd "cactus-bin-v2.0.4.tar.gz"
source venv/bin/activate
export PATH=$(pwd)/bin:$PATH
export PYTHONPATH=$(pwd)/lib:$PYTHONPATH
```
Можно создать файл activate.sh в папке venv, чтобы не писать три команды отдельно:

```bash
vim activate.sh
```
 Внутри файла пишем:
 ```bash
 source venv/bin/activate
 export PATH=$(pwd)/bin:$PATH
 export PYTHONPATH=$(pwd)/lib:$PYTHONPATH
 ```
 Далее активируем окружение каткуса с уже указанными путями через команду:
 ```bash
 source venv/activate.sh
 ```

Файл конфига:

```text
((simHuman_chr6:0.144018,(simMouse_chr6:0.084509,simRat_chr6:0.091589)mr:0.271974):0.020593,(simCow_chr6:0.18908,simDog_chr6:0.16303):0.032898);

simCow_chr6 https://raw.githubusercontent.com/UCSantaCruzComputationalGenomicsLab/cactusTestData/master/evolver/mammals/loci1/simCow.chr6
simDog_chr6 https://raw.githubusercontent.com/UCSantaCruzComputationalGenomicsLab/cactusTestData/master/evolver/mammals/loci1/simDog.chr6
simHuman_chr6 https://raw.githubusercontent.com/UCSantaCruzComputationalGenomicsLab/cactusTestData/master/evolver/mammals/loci1/simHuman.chr6
simMouse_chr6 https://raw.githubusercontent.com/UCSantaCruzComputationalGenomicsLab/cactusTestData/master/evolver/mammals/loci1/simMouse.chr6
simRat_chr6 https://raw.githubusercontent.com/UCSantaCruzComputationalGenomicsLab/cactusTestData/master/evolver/mammals/loci1/simRat.chr6
```

Это дерево первой строкой.

Потому пробел.

Потом строки название сборки пробел путь к фасте.

Не надо пытаться давать параметры ядер, это ломает:

```
--defaultCores 192 --defaultMemory 1000G
```
