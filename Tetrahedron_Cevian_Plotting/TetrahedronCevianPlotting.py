
import os
from collections import defaultdict

import numpy as np

from tqdm import tqdm

import plotly as py
import plotly.graph_objects as go

### CONSTANTS ###

# verticies of equilateral triangle along unit 1 circle
SQRT3_OVER2_ = np.sqrt(3, dtype=np.float64)/2.
EQUILATERALVERTICIES = np.array([np.array([.0, 1.0], dtype=np.float64),
                                 np.array([-1 * SQRT3_OVER2_, -1/2],
                                          dtype=np.float64),
                                 np.array([1 * SQRT3_OVER2_, -1/2],
                                          dtype=np.float64)
                                 ]
                                )

# verticies of tetrahedron along unit 1 sphere
SQRT_2OVER3_ = np.sqrt(2/3, dtype=np.float64)
SQRT_8OVER9_ = np.sqrt(8/9, dtype=np.float64)
SQRT_2OVER9_ = np.sqrt(2/9, dtype=np.float64)

TETRAHEDRONVERTICIES = np.array([np.array([0, 0, 1],
                                          dtype=np.float64),
                                 np.array([SQRT_8OVER9_, 0, -1/3],
                                          dtype=np.float64),
                                 np.array(
                                     [-1*SQRT_2OVER9_, SQRT_2OVER3_, -1/3], dtype=np.float64),
                                 np.array(
                                     [-1*SQRT_2OVER9_, -1*SQRT_2OVER3_, -1/3], dtype=np.float64)
                                 ])


class TetrahedronCevianPlotting:
    '''General transformation and plotting class for phase diagrams.

    Phase plots display the change in the states of a dynamic system with
    respect to themselves.
    for example the state of a system with states:
        [s0, s1, s2, ..... sx] with n = timestep number
    could be ploted in 2 dimensions by iteratively vaying the state values
    along the x & y axes in the following manner:
                   x ,    y
        p0      = [s0,    s1]
        p1      = [s1,    s2]
        p2      = [s2,    s3]
                .
                .
                .
        p(n-1)  = [s(n-1), s(n)]

    The same can be applied to the same data but in 2 dimensions by increasing
    the axes count to 3 (x, y, z) and applying the same technique above:

                   x ,    y
        p0      = [s0,    s1,     s2]
        p1      = [s1,    s2,     s3]
        p2      = [s2,    s3,     s4]
                .
                .
                .
        p(n-2)  = [s(n-2), s(n-1), s(n)]

        note how the number of points decreases given the # of axes used

    This however limits the number of axes available for plotting given the
    maximum of 3 dimensions available for visualisation with standard
    cartesian coordinates.

    A re-interpretation of this process in triangular space allows for plotting
    of higher dimensional data in the 3 cartesian dimensions available.

    In the case of the 3D phase ploting, the data is in the form of:
        p(n-2)  = [s(n-2), s(n-1), s(n)]

    If the 3 data points are normalised by finding the weight of every state
        w(n-2) = [s(n-2), s(n-1), s(n)] / sum(p(n-2))

        it is evident that 
        sum(w(n-2)) = 1

    As such, the values of w(n-2) can be interpreted as the proportional
    weights of p(n-2) that can represent the cevian proportions of a point
    plotted inside the equilateral triangle given that the weighted cevian 
    coordinates of the verticies of a triangle correspond to:
        (1,0,0), (0,1,0), (0,0,1) with the centroid at (1/3, 1/3, 1/3)

    Since these three points are representable on the standard (x,y) plane,
    for example by considering them as the points of an equilateral triangle
    laying on a unit circle, the reinterpreted data is effectively mapped from a
    3D space onto a 2D one.

    The same holds for the next higher order equivalent. In this case a 4D set
    of phase points, can be represented as weighted points within a tetrahedron
    and as such are visualisable in standard 3D cartesian coordiantes. In this
    case a the circumcribed sphere of a tetrahedron is used with cevian
    verticies:
        (1,0,0,0), (0,1,0,0), (0,0,1,0), (0,0,0,1)
        with the centroid at (1/4, 1/4, 1/4, 1/4)

    This class transfroms an arbitrary set of generated state values for an 
    arbitrary number of states, and plots each of these value sets onto either
    a 2D representation of 3 phase state length (in a triangle) or a 3D
    representation of 4 phase state length (in a tetrahedron)

    Attributes
    ----------
    pdim: int
        the required cevian dimension to transform to
    ref_v : list
        the verticies used for plotting in cartesian coords
    outfile : str
        plotting output file absolute path
    overlay_plots : bool
        overlay indpendent state data or generate separate plots for each

    Methods
    -------
    get_plotting_coords(raw_vals)
        transforms appropriate single set of appropriate length phase states
        data into cartesian coords
    process_data(pdata)
        process an entire data set of arbitrary length state values for an
        arbitrary number of different states
    plot_data(pdata)
        use process_data to transform teh data and generate teh required plot(s)
    '''

    allowed_ptypes = {'triangle': 3,
                      'tetrahedron': 4
                      }

    verticies_ref = {'triangle': EQUILATERALVERTICIES,
                     'tetrahedron': TETRAHEDRONVERTICIES
                     }

    def __init__(self, ptype, fname, fdir=None, overlay_plots=True):
        '''
        Parameters
        ----------
        ptype : str
            Type of plot requested (triangle or tetrahedron)
        fname : str
            Properly formatted file name
        fdir : str, optional
            Properly formatted output directory
            note: if not passed, default directory used
        overlay_plots: bool
            Overlay the various states data onto a single plot or generate
            separate plots for each
        '''
        self.ptype = ptype
        self.pdim = ptype
        self.ref_v = ptype
        self.outfile = [fdir, fname]
        self.overlay_plots = overlay_plots

    @property
    def pdim(self):
        return self.__pdim

    @pdim.setter
    def pdim(self, ptype):
        try:
            self.__pdim = self.allowed_ptypes[ptype]
        except:
            raise ValueError(f'ptype must be one of \
                               {self.allowed_ptypes.keys()}')

    @property
    def ref_v(self):
        return self.__ref_v

    @ref_v.setter
    def ref_v(self, ptype):
        self.__ref_v = self.verticies_ref[ptype]

    @property
    def outfile(self):
        return self.__outfile

    @outfile.setter
    def outfile(self, fdata):
        fdir, fname = fdata
        # check if dir exists, if not recursively make it
        if not fdir:
            fdir = os.getcwd() + '\\visualisations'
            if not os.path.exists(fdir):
                os.makedirs(fdir)
        elif not os.path.isabs(fdir):
            fdir = os.path.dirname(os.path.abspath('\\' + fdir))

        self.__outfile = fdir
        self.__outfile += '\\' + fname + '.html'

    def get_plotting_coords(self, rvals):
        '''Generate the cartesian plotting coordinates for a set of state values
        based on the chosen plotting type (triangle, tetrahedron)

        The raw_values must be of the appropriate length required for the 
        chosen plotting type (e.g. 3 state values to generate 2D trangle plot
        & 4 state values required to generate 3D tetrahedron plot)

        Parameters
        ----------
        rvals : np.ndarray
            raw state value to be processed. shape of array must (n,) where
            n = appropriate number of values

        Raises
        ------

        '''

        cart_point = None
        for idx in range(self.pdim-1):
            try:
                prop_weight = rvals[1]/(rvals[0] + rvals[1])
            except:
                prop_weight = 0.0

            # initial run ?
            if not isinstance(cart_point, np.ndarray):
                to_vec = self.ref_v[1] - self.ref_v[0]
                cart_point = self.ref_v[0] + prop_weight * to_vec
            else:
                to_vec = self.ref_v[idx + 1] - cart_point
                cart_point = cart_point + prop_weight * to_vec

            rvals = np.concatenate((np.sum(rvals[:2]), rvals[2:]),
                                   axis=None
                                   )
        return cart_point

    def process_data(self, rdata):
        '''Process an entire set of arbitrary length state values for an 
        arbitrary number of states

        Parameters
        ----------
        rdata : np.ndarray
            raw state values data to be processed. shape of array must (n,m)
            where n = count of state values & m = number of independent states

        Raises
        ------

        '''

        # cartesian point values
        p_values = defaultdict()

        # check for number of independent states in the data
        # each state is the equivalent to an independent plot trace
        try:
            points_count, traces_count = rdata.shape
        except:
            rdata = np.array([rdata]).T
            points_count, traces_count = rdata.shape

        # number of return values
        itts_count = points_count - self.pdim + 1
        print('Processing Traces:')
        for trace_num in tqdm(range(traces_count)):
            # stack data for processing
            # i.e. generate [s0, s1 ... s(pdim)] arrays
            stacked_values = tuple(rdata[idx:idx+itts_count, trace_num]
                                   for idx in range(self.pdim)
                                   )
            stacked_values = np.column_stack(stacked_values)

            # apply conversion to get cartesian points
            cartesian_points = np.apply_along_axis(func1d=self.get_plotting_coords,
                                                   axis=1,
                                                   arr=stacked_values
                                                   )

            cartesian_points = np.unique(cartesian_points, axis=0)

            # map values onto independent cartesian plotting axes
            axes = 'xyz'
            p_values[trace_num] = {axes[idx]: cartesian_points[:, idx]
                                   for idx in range(self.pdim-1)
                                   }
        return dict(p_values)

    def get_colorscale(self):
        colorscale = ['#000086',
                      '#000086',
                      '#000086',
                      '#000086',
                      '#000086',
                      '#000086',
                      '#33198d',
                      '#4f2f94',
                      '#674398',
                      '#b23978',
                      '#c95366',
                      '#de6c47',
                      '#f28500',
                      '#918a42',
                      '#a59659',
                      '#b4a36f',
                      '#c2b186',
                      '#cdc09e',
                      '#d7cfb5',
                      '#dfdfcd',
                      '#000086',
                      '#33198d',
                      '#4f2f94',
                      '#674398',
                      '#b23978',
                      '#c95366',
                      '#de6c47',
                      '#f28500',
                      '#918a42',
                      '#a59659',
                      '#b4a36f',
                      '#c2b186',
                      '#cdc09e',
                      '#d7cfb5',
                      '#dfdfcd',
                      ]

        colorscale = [[val, color]
                      for val, color in zip(np.linspace(start=0.0,
                                                        stop=1.0,
                                                        num=len(colorscale)),
                                            colorscale
                                            )
                      ]

        return colorscale

    def generate_plots(self, pdata, tlabels):
        '''Generate required plots for provided processed/mapped data

        Parameters
        ----------
        pdata : dict
            processed (weighted & mapped) state values data to be plotted.
            dict structure: 
                pdata[trace_number][axis_title: (x,y) for 2D & (x,y,z) for 3D]


        Raises
        ------

        '''

        plot_traces = []
        fig = None
        colorscale = self.get_colorscale()
        print('Generating Traces: ')

        if self.ptype == 'triangle':
            x_vals, y_vals, colors, text = [], [], [], []

            for t_id, data in tqdm(pdata.items()):
                tcolor = tlabels[t_id]

                for x_val, y_val in zip(pdata[t_id]['x'], pdata[t_id]['y']):
                    x_vals.append(x_val)
                    y_vals.append(y_val)
                    colors.append(tcolor)

            plot_traces.append(go.Scattergl({'name': tlabels[t_id],
                                             'mode': 'markers',
                                             'x': x_vals,
                                             'y': y_vals,
                                             'marker': {'symbol': 'circle',
                                                        'color': colors,
                                                        'colorbar': {'title': 'r-value'},
                                                        'colorscale': colorscale,
                                                        'size': 2.5,
                                                        'opacity': 0.75
                                                        }
                                             }
                                            )
                               )

            v_x, v_y, v_text = [], [], []
            v_coords = [(1, 0, 0), (0, 1, 0), (0, 0, 1)]
            for v_idx, vert in enumerate(self.ref_v):
                v_x.append(vert[0])
                v_y.append(vert[1])
                v_text.append(f'X @ t = {v_idx}<br>{v_coords[v_idx]}')

            plot_traces.append(go.Scattergl({'mode': 'lines+markers+text',
                                             'x': v_x,
                                             'y': v_y,
                                             'text': v_text,
                                             'marker': {'symbol': 'circle',
                                                        'color': 'rgb(245, 245, 245)',
                                                        'size': 8,
                                                        'opacity': 1.0
                                                        }
                                             }
                                            )
                               )
            fig = go.Figure(data=plot_traces)

            fig.update_layout({'autosize': False,
                               'showlegend': False,
                               'width': 1600,
                               'height': 1600,
                               'margin': go.layout.Margin(l=20,
                                                          r=20,
                                                          b=0,
                                                          t=0,
                                                          pad=0
                                                          ),
                               'paper_bgcolor': 'rgba(64,64,64,0.5)',
                               'plot_bgcolor': 'rgba(64,64,64,0.5)',
                               'scene': {'xaxis': {'showbackground': False},
                                         'yaxis': {'showbackground': False}
                                         }
                               }
                              )

        else:
            x_vals, y_vals, z_vals, colors, text = [], [], [], [], []

            for t_id, data in tqdm(pdata.items()):
                tcolor = tlabels[t_id]

                for x_val, y_val, z_val in zip(pdata[t_id]['x'],
                                               pdata[t_id]['y'],
                                               pdata[t_id]['z']):
                    x_vals.append(x_val)
                    y_vals.append(y_val)
                    z_vals.append(z_val)
                    colors.append(tcolor)

            plot_traces.append(go.Scatter3d({'mode': 'markers',
                                             'x': x_vals,
                                             'y': y_vals,
                                             'z': z_vals,
                                             'marker': {'symbol': 'circle',
                                                        'color': colors,
                                                        'colorbar': {'title': 'r-value'},
                                                        'colorscale': colorscale,
                                                        'size': 2.5,
                                                        'opacity': 0.75
                                                        }
                                             }
                                            )
                               )

            # add vericies traces
            v_x, v_y, v_z, v_text = [], [], [], []
            v_coords = [(1, 0, 0, 0), (0, 1, 0, 0), (0, 0, 1, 0), (0, 0, 0, 1)]
            for v_idx, vert in enumerate(self.ref_v):
                v_x.append(vert[0])
                v_y.append(vert[1])
                v_z.append(vert[2])
                v_text.append(f'X @ t = {v_idx}<br>{v_coords[v_idx]}')

            plot_traces.append(go.Scatter3d({'mode': 'lines+markers+text',
                                             'x': v_x,
                                             'y': v_y,
                                             'z': v_z,
                                             'text': v_text,
                                             'marker': {'symbol': 'circle',
                                                        'color': 'rgb(245, 245, 245)',
                                                        'size': 10,
                                                        'opacity': 1.0
                                                        }
                                             }
                                            )
                               )

            fig = go.Figure(data=plot_traces)

            fig.update_layout({'autosize': False,
                               'showlegend': False,
                               'width': 1600,
                               'height': 1600,
                               'margin': go.layout.Margin(l=20,
                                                          r=20,
                                                          b=0,
                                                          t=0,
                                                          pad=0
                                                          ),
                               'paper_bgcolor': 'rgba(16,16,16,0.75)',
                               'plot_bgcolor': 'rgba(16,16,16,0.75)',
                               'scene': {'xaxis': {'showbackground': False},
                                         'yaxis': {'showbackground': False},
                                         'zaxis': {'showbackground': False}
                                         }
                               }
                              )
        print(f'Writing Plot to disk at: {self.outfile}')
        py.offline.plot(fig, filename=self.outfile, auto_open=False)

    def plot_data(self, rawdata, labels):
        '''Generate required plots for the given data

        Parameters
        ----------
        rawdata : np.ndarray
            raw state values data to be processed. shape of array must (n,m)
            where n = count of state values & m = number of independent states

        Raises
        ------

        '''

        # Process raw data to get cartesian points
        plotting_values = self.process_data(rawdata)
        self.generate_plots(plotting_values, labels)
