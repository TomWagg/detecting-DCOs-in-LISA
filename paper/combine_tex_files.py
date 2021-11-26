with open("combined.tex", "w") as combined:
    with open("tex/paper.tex", "r") as paper:
        for line in paper:
            if line.strip()[:7] == "\\input{":
                with open("tex/" + line.strip()[7:].replace("}", "").replace("\n", ""), "r") as section:
                    combined.writelines(section.readlines())
            else:
                combined.write(line.replace("../figures/", "./figures/"))
