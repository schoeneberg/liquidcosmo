import matplotlib as mpl
def matplotlib_defaults(rwth_colors=3, backend=None, update_tex=True, update_size = [10,12,14]):

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

  # Set colors
  rwth={}
  RWTHcolors={}
  RWTHcolors_DCH = {}

  # RWTH CD colors
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
  rwth_array = [rwth[c] for c in ['Blau','Grun','Rot','Hellblau','Lila','Schwarz','Magenta','Gelb','Petrol','Maigrun','Orange','Bordeaux','Turkis','Violett']]

  RWTH_goods = {'Red':'#E37C80','Blue':'#7A98F6','Green':'#88B27A','Orange':'#F3BE82','Grey':'#ABABAB','Purple':'#B87294'}
  rwth_array2 = [RWTH_goods[c] for c in ['Red','Blue','Green','Orange','Purple','Grey']]

  RWTHcolors = {'Red':['#E37C80','#CE121F'],'Blue':['#7A98F6','#1157EF'],'Green':['#88B27A','#297C09'],'Orange':['#F3BE82','#ED920F'],'Grey':['#ABABAB','#737373'],'Purple':['#B87294','#88004C']}
  RWTHcolors_DCH = {'Red':['#E69679','#CC071E'],'Blue':['#8EBAE5','#00549F'],'Green':['#B8D698', '#57AB27'],'Orange':['#FDD48F','#F6A800'], 'Purple':['#BCB5D7','#7A6FAC'],'Grey':['#ABABAB','#737373'],'Black':['#646567','#9C9E9F'],'Turk':['#0098A1','#89CCCF']}  
  rwth_array3 = [RWTHcolors_DCH[c][0] for c in ['Red','Blue','Green','Orange','Purple','Grey','Turk','Black']]
  if rwth_colors == 1:
    mpl.rcParams['axes.prop_cycle'] = mpl.cycler(color=rwth_array) 
  elif rwth_colors == 2:
    mpl.rcParams['axes.prop_cycle'] = mpl.cycler(color=rwth_array2)
  elif rwth_colors == 3:
    mpl.rcParams['axes.prop_cycle'] = mpl.cycler(color=rwth_array3)
  else:
    raise Exception("rwth_colors option can only be 1 or 2 or 3")

