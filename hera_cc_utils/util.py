""" Random stuff. """

@ticker.FuncFormatter
def degree_ticks(x, pos):
    return r"%.0f$^\circ$" % x

def top_cbar(ax, cax, label='Jy/beam', size='4%', pad=0.1, length=5,
    labelsize=16, fontsize=20, minpad=1):

    divider = make_axes_locatable(ax)
    cbax = divider.append_axes('top', size=size, pad=pad)
    cbar = fig.colorbar(cax, cax=cbax, orientation='horizontal')
    cbax.grid(False)
    cbax.xaxis.set_ticks_position('top')
    cbax.tick_params(length=length, labelsize=labelsize)
    cbax.xaxis.set_label_position('top')
    cbax.set_xlabel(label, fontsize=fontsize)
    return cbax, cbar
