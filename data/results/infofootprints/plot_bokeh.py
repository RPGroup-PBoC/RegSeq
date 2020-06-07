import bokeh
from bokeh.io import output_file, show
from bokeh.models.widgets import Dropdown
from bokeh.models.widgets import Select
import pandas as pd
import numpy as np
from bokeh.layouts import layout
from bokeh.plotting import figure, show
bokeh.plotting.output_file("/home/bill/panel-examples/test.html")
p = figure(x_range=[str(x) for x in range(160)], width=200, height=30)
p.vbar(x=[str(x) for x in range(160)], top=data_0,width=.9)

z = layout([
  [select, p],
], sizing_mode='stretch_both')


from bokeh.models import CustomJS, Slider, ColumnDataSource
from bokeh.plotting import figure, output_file, show

output_file("dropdown.html")


select = Select(title="Option:", value="foo", options=["foo", "bar", "baz", "quux"])

data = pd.io.parsers.read_csv('all_info')
data_1 = data.loc[0,[str(x) for x in range(160)]]
print(data_1)

p = figure(x_range=[str(x) for x in range(160)], width=200, height=30)
p.vbar(x=[str(x) for x in range(160)], top=data_1,width=.9)

z = layout([
  [select, p],
], sizing_mode='stretch_both')

bokeh.io.save(z)

