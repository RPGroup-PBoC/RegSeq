# +
import matplotlib.pyplot as plt 
import matplotlib as mpl
import seaborn as sns
import matplotlib
from matplotlib import rcParams

rcParams['axes.titlepad'] = 20 

# Default RP plotting style

# Default RP plotting style
def pboc_style_mpl():
    """
    Formats matplotlib plotting enviroment to that used in 
    Physical Biology of the Cell, 2nd edition.
    """
    
    tw = 1.5
    rc = {'lines.linewidth': 1.25,
          'axes.labelsize': 14,
          'axes.titlesize': 18,
          'axes.facecolor': '#E3DCD0',
          
          'xtick.major' : 16,
          'ytick.major' : 16,
          'xtick.major.width': tw,
          'xtick.minor.width': tw,
          'ytick.major.width': tw,
          'ytick.minor.width': tw,
          'xtick.labelsize': 'large',
          'ytick.labelsize': 'large',
          'xtick.labelsize': 10,
          'ytick.labelsize': 10,
          'font.family': 'sans-serif',
          'font.sans-serif': 'Lucida Grande',
          'grid.linestyle': '-',
          'grid.linewidth': 0.5,
          'grid.color': '#ffffff',
          'legend.fontsize': 10,
          'figure.dpi': 300,
          'savefig.dpi': 300}
    plt.rc('text.latex', preamble=r'\usepackage{sfmath}')
    plt.rc('xtick.major', pad=-1)
    plt.rc('ytick.major', pad=-1)
    plt.rc('mathtext', fontset='stixsans', sf='sansserif')
    plt.rc('figure', figsize=[3.5, 2.5])
    plt.rc('svg', fonttype='none')
    plt.rc('legend', title_fontsize='8', frameon=True, 
           facecolor='#E3DCD0', framealpha=1)
    sns.set_style('darkgrid', rc=rc)
    sns.set_palette("colorblind", color_codes=True)
    sns.set_context('notebook', rc=rc)


def pboc_style_mpl_2():
      
    tw = 1.5

    rc = {'lines.linewidth': 2,
        'axes.labelsize': 18,
        'axes.titlesize': 21,
        'xtick.major' : 16,
        'ytick.major' : 16,
        'xtick.major.width': tw,
        'xtick.minor.width': tw,
        'ytick.major.width': tw,
        'ytick.minor.width': tw,
        'xtick.labelsize': 'large',
        'ytick.labelsize': 'large',
        'font.family': 'sans',
        'weight':'bold',
        'grid.linestyle': ':',
        'grid.linewidth': 1.5,
        'grid.color': '#ffffff',
        'mathtext.fontset': 'stixsans',
        'mathtext.sf': 'fantasy',
        'legend.frameon': True,
        'legend.fontsize': 12, 
       "xtick.direction": "in","ytick.direction": "in"}



    plt.rc('text.latex', preamble=r'\usepackage{sfmath}')
    plt.rc('mathtext', fontset='stixsans', sf='sans')
    sns.set_style('ticks', rc=rc)

    #sns.set_palette("colorblind", color_codes=True)
    sns.set_context('notebook', rc=rc)


def pboc_style_bokeh():
    '''
    Formats bokeh plotting enviroment to that used in 
    Physical Biology of the Cell, 2nd edition.
    '''
    theme_json = {'attrs':{'Axis': {
            'axis_label_text_font': 'Helvetica',
            'axis_label_text_font_style': 'normal'
            },
            'Legend': {
                'border_line_width': 1.5,
                'background_fill_alpha': 0.5
            },
            'Text': {
                'text_font_style': 'normal',
               'text_font': 'Helvetica'
            },
            'Title': {
                #'background_fill_color': '#FFEDC0',
                'text_font_style': 'normal',
                'align': 'center',
                'text_font': 'Helvetica',
                'offset': 2,
            }}}

    return theme_json
# -


