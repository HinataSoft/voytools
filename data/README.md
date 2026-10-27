# data/ – transliterace použité v repu

Kopie souborů ze stránky Reného Zandbergena <https://www.voynich.nu/transcr.html> (stav 05/2025),
formát IVTFF 2.0 (LSI: 1.5). Jsou tu proto, aby cloudové session nepotřebovaly síť mimo GitHub.
Soubory jsou dílem jejich autorů (Zandbergen & Landini, Claston, Takahashi, Stolfi a dal.) a jsou
zde jen pro výzkumné použití; při případném zveřejnění repa zvažte, zda je nahradit skriptem
`download.sh` (stáhne je znovu z voynich.nu, pokud to síť prostředí dovolí).

| Soubor | Co to je | Abeceda |
|---|---|---|
| `ZL3b-n.txt` | Zandbergen–Landini v. 3b – **primární pracovní soubor** | Eva |
| `GC2a-n.txt` | Glen Claston v. 2a – nezávislé úplné čtení, **kontrola** | v101 |
| `RF1b-er.txt` | Referenční (kombinace ZL+GC), základní Eva | Eva |
| `IT2a-n.txt` | Takahashi (stav 1999) | Eva-T |
| `LSI_ivtff_0d.txt` | Landini–Stolfi interlineár (více transkribentů, kód za `;`) | Eva-T |

Kontrolní součty (počet stránek / lokusů), jak je počítá `python -m voytools summary`:
ZL3b 227 / 5 385, GC2a 226 / 5 367.

Úplný přehled všech transliterací, abeced a formátu je v `D:\Work\Voynich\transliterations\README.md`
(mimo repo). Definice formátu: <https://www.voynich.nu/software/ivtt/IVTFF_format.pdf>.
