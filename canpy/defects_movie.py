# -*- coding: utf-8 -*-
"""
Created on Sun Feb 14 20:43:12 2016

@author: Adrien
"""

import canpy.point_defects


def defects_movie(dl,t_start='standard',t_end='standard'):
    '''
    Produce a movie in mp4 format of the defects.
    t_start is the time at which we want the movie to begin and t_end is the time
    at which we want the movie to end.
    If t_start and t_end are kept at the 'standard' value the movie will start
    at 0ps and end at the last available snapshot. The user is free to change the
    beginning and finishing time. If so these times have to be provided in ps.
    '''
    
    if t_start == 'standard':
        snap_start=0
    else:
        snap_start = round(t_start/(dl.inter*dl.time_step))
    
    if t_end == 'standard':
        snap_end=dl.nb_snapshots-1
    else:
        snap_end = round(t_end/(dl.inter*dl.time_step))
    
    i=snap_start
    while i <= snap_end:
        a = package.point_defects(dl,i*dl.inter*dl.time_step)
        
        
        
        
        
        
    fig = plt.figure()
    ax = fig.gca(projection='3d')
        
        
        
    def init():
        
        ax.set_xlabel('x axis[Angström]')
        ax.set_ylabel('y axis[Angström]')
        ax.set_zlabel('z axis[Angström]')
            
        ax.set_xlim3d(0,self.length)
        ax.set_ylim3d(0,self.length)
        ax.set_zlim3d(0,self.length)
            
        ax.legend(frameon = True, fancybox = True, ncol = 1, fontsize = 'x-small', loc = 'lower right')
            
        
        
    def animate(i):
        
        
        
        ax.scatter(xi, yi, zi, label='Interstitials', c='r', marker='^')
        ax.scatter(xv, yv, zv, label='Vacancies', c='b', marker='o')

        ax.set_title('Interstitials and vacancies at '+str(self.time)+' ps')
        
    # Animate
    anim = animation.FuncAnimation(fig, animate, init_func=init,frames=360, interval=speed*20, blit=True)
    # Save
    anim.save('rot_frame_anim_'+str(self.time)+'ps.mp4', fps=30, extra_args=['-vcodec', 'libx264'])
    