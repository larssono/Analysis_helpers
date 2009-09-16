import urllib
from PIL import Image
import StringIO

def venn_diagram_url(data, width=400, height=200):
    a,b,c = data.values()
    values = [width, height, len(a), len(b), len(c)]
    values.append(len(a & b))
    values.append(len(a & c))
    values.append(len(c & b))
    values.append(len(a & b & c))
    values.append('|'.join(data.keys()))
    base_str = 'http://chart.apis.google.com/chart?'
    args_str = 'cht=v&chs=%sx%s&chd=t:%s,%s,%s,%s,%s,%s,%s&chdl=%s'   
    return base_str + (args_str % tuple(values))

data={'anova': set((1, 2, 3,4,5)),
      'manova': set((2, 4)),
      'hosvd': set((3, 4, 5))}

img = StringIO.StringIO(urllib.urlopen(venn_diagram_url(data)).read())
Image.open(img).show()
