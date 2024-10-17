import matplotlib as mpl

default_colors = None
plot_settings = None

def get_colors():
  colors = {}
  rwth = {}
  rwth['Blau']=(0./255,84./255,159./255)
  rwth['Hellblau']=(142./255,186./255,229./255)
  rwth['Schwarz']=(0./255,0./255,0./255)
  rwth['Magenta']=(227./255,0./255,102./255)
  rwth['Gelb']=(255./255,237./255,0./255)
  rwth['Petrol']=(0./255,97./255,101./255)
  rwth['Turkis']=(0./255,152./255,161./255)
  rwth['Grun']=(87./255,171./255,39./255)
  rwth['Maigrun']=(189./255,205./255,0./255)
  rwth['Orange']=(246./255,168./255,0./255)
  rwth['Rot']=(204./255,7./255,30./255)
  rwth['Bordeaux']=(161./255,16./255,53./255)
  rwth['Violett']=(97./255,33./255,88./255)
  rwth['Lila']=(122./255,111./255,17./255)
  colors['rwth'] = rwth

  goods = {'Red':'#E37C80','Blue':'#7A98F6','Green':'#88B27A','Orange':'#F3BE82','Grey':'#ABABAB','Purple':'#B87294'}
  colors['v0'] = goods

  MPcolors = {'Red':['#E37C80','#CE121F'],'Blue':['#7A98F6','#1157EF'],'Green':['#88B27A','#297C09'],'Orange':['#F3BE82','#ED920F'],'Grey':['#ABABAB','#737373'],'Purple':['#B87294','#88004C']}
  colors['MP'] = MPcolors

  # Colors from Deanna C Hooper
  DCH_colors = {'Red':['#E69679','#CC071E'],'Blue':['#8EBAE5','#00549F'],'Green':['#B8D698', '#57AB27'],'Orange':['#FDD48F','#F6A800'], 'Purple':['#BCB5D7','#7A6FAC'],'Grey':['#ABABAB','#737373'],'Black':['#646567','#9C9E9F'],'Turk':['#0098A1','#89CCCF']}
  colors['DCH'] = DCH_colors

  # Works quite well in general, has a nice look to it. Also looks nice in greyscale, actually.
  v1_colors = {'Red':'#C9031A','Blue':'#00519E','Green':'#57AB26', 'Orange':'#F5A600', 'Purple':'#786EAB', 'Grey':'#919191', 'DarkBlue':'#000688'}
  colors['v1'] = v1_colors

  # Specifically for colorblind people
  v1_CB_colors = {'Blue':'#377eb8', 'Orange':'#ff7f00', 'Green':'#4daf4a', 'Pink':'#f781bf', 'Brown':'#a65628', 'Purple':'#984ea3','Grey':'#999999', 'Red':'#e41a1c', 'Yellow':'#dede00'}
  colors['v1_CB'] = v1_CB_colors

  return colors

# Make available to import (!)
colors = get_colors()

def matplotlib_defaults(selection='v1', backend=None, update_tex=True, update_size = [10,12,14]):
  global default_colors

  # If required, set backend
  if backend:
    mpl.use(backend)

  if update_tex:
    mpl.rc('font',**{'family':'serif','serif':['Computer Modern Roman']})
    mpl.rc('text', usetex=True)
    mpl.rc('text.latex', preamble=r'\usepackage{amsmath}')
    mpl.rc('image', cmap='viridis')

  if update_size:
    try:
      success = (len(update_size) == 3)
    except:
      success = False
    if success:
      SMALL_SIZE = update_size[0]
      MEDIUM_SIZE = update_size[1]
      LARGE_SIZE = update_size[2]
      mpl.rc('font', size=MEDIUM_SIZE)          # controls default text sizes
      mpl.rc('axes', titlesize=MEDIUM_SIZE)     # fontsize of the axes title
      mpl.rc('axes', labelsize=LARGE_SIZE)    # fontsize of the x and y labels
      mpl.rc('xtick', labelsize=SMALL_SIZE)    # fontsize of the tick labels
      mpl.rc('ytick', labelsize=SMALL_SIZE)    # fontsize of the tick labels
      mpl.rc('legend', fontsize=MEDIUM_SIZE)    # legend fontsize
      mpl.rc('figure', titlesize=LARGE_SIZE)  # fontsize of the figure title
    else:
      print("UNABLE TO SET SIZES FOR MATPLOTLIB")

  colors = get_colors()
  
  choice = {1:'v1',2:'v1_CB',3:'MP',4:'DCH',5:'rwth'}
  if isinstance(selection,int) and selection>0 and selection<len(choice)+1:
    cname = choice[selection]
  elif selection in choice.values():
    cname = selection
  else:
    raise ValueError("Unknown color choice '{}'".format(selection))

  if cname == 'rwth':
    carray = {c:colors[cname][c] for c in ['Blau','Grun','Rot','Hellblau','Lila','Schwarz','Magenta','Gelb','Petrol','Maigrun','Orange','Bordeaux','Turkis','Violett']}
  elif cname == 'MP':
    carray = {c:colors[cname][c] for c in ['Red','Blue','Green','Orange','Purple','Grey']}
  elif cname=='DCH':
    carray = {c:colors[cname][c][0] for c in ['Red','Blue','Green','Orange','Purple','Grey','Turk','Black']}
  else:
    carray = colors[cname]

  selected_colors = list(carray.values())

  mpl.rcParams['axes.prop_cycle']  = mpl.cycler(color=selected_colors)
  default_colors = colors[cname]

  return selected_colors

def initialize_plots(legend_frame = True, legend_fontsize = 15, axes_fontsize = 15, axes_labelsize = 20, **kwargs):
  global plot_settings
  my_colors = matplotlib_defaults(**kwargs)
  import getdist.plots
  gdplotsettings = getdist.plots.GetDistPlotSettings()
  gdplotsettings.solid_colors = my_colors
  gdplotsettings.line_styles = [("-",c) for c in my_colors]
  gdplotsettings.solid_contour_palefactor = 0.5
  gdplotsettings.linewidth_contour = 2.0
  gdplotsettings.linewidth = 2.0
  gdplotsettings.legend_fontsize = legend_fontsize
  gdplotsettings.axes_fontsize = axes_fontsize
  gdplotsettings.axes_labelsize = axes_labelsize
  plot_settings = gdplotsettings

  if not legend_frame:
    g.settings.figure_legend_frame = False

  return gdplotsettings

# Little trick to keep the other modules always referencing this exact instance of the plot settings
def get_plot_settings():
  global plot_settings
  return plot_settings
