

def parse_table(p, table_tag):
    pytable = []
    for tr in table_tag('tr'):
        pytable.append([])

        for th in p(tr)('th'):
            pytable[-1].append(p(th).text())

        for td in p(tr)('td'):
            pytable[-1].append(p(td).text())
    return pytable
