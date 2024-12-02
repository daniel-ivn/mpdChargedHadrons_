## Установить pyROOT и запустить postprocess

```
PYTHONPATH=/home/dasha/NICA/mpdGit/venv/bin/python
PYTHONPATH=$PYTHONPATH:$ROOTSYS/lib
source thisroot.sh (root6)
pip install rootpy
python postprocess.py -i 'input/mpdpid10.root' -d output -s input/settings.json
```

## Файлы

 __input/train_output/\*__ - выходные файлы паровоза Киреева, взаты из /junk/kireev/paper_mpdpid*.root

**input/postprocess_mpdpid10.root** - постпроцесс файлов паровоза Киреева из train_output

**input/nuclei_spectra** - то, что я сама запускала раньше на маленькой статистике (на 1 файле.)

**ChemicalPotential** - строит фазовую диаграмму и отношения





