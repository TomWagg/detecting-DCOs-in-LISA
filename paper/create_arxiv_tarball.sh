rm -rf arxiv/
mkdir arxiv
mkdir arxiv/figures
cp figures/fig*.pdf arxiv/figures
cp figures/fig*.png arxiv/figures
cp ./{tex/*.bst,tex/*.cls,tex/ORCID-ID.png,tex/*.bib} arxiv
cp bbl.txt arxiv/paper.bbl
python combine_tex_files.py
mv combined.tex arxiv/paper.tex
tar -czvf arxiv.tar.gz arxiv/
rm -rf arxiv/
