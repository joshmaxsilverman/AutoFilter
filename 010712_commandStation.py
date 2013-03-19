import os
import pprint
import random
import wx
import csv
import numpy

# 1 - The recommended way to use wx with mpl is with the WXAgg backend. 

import matplotlib
matplotlib.use('WXAgg')
from matplotlib.figure import Figure
from matplotlib.backends.backend_wxagg import \
    FigureCanvasWxAgg as FigCanvas, \
    NavigationToolbar2WxAgg as NavigationToolbar
from matplotlib import pyplot as plt

class MainFrame(wx.Frame):
	
	title = "Function View"
	save_csv_option = 0
	
	def __init__(self):
		wx.Frame.__init__(self, None, -1, self.title)
		
		self.outputName = "file.txt"
		
		self.create_menu()
		self.create_status_bar()
		self.create_main_panel()
		
		self.draw_figure()
		
	def create_menu(self):
	
		self.menubar = wx.MenuBar()
        
        menu_file = wx.Menu()
        m_expt = menu_file.Append(-1, "&Save plot\tCtrl-S", "Save plot to file")
        self.Bind(wx.EVT_MENU, self.on_save_plot, m_expt)
        menu_file.AppendSeparator()
        m_exit = menu_file.Append(-1, "E&xit\tCtrl-X", "Exit")
        self.Bind(wx.EVT_MENU, self.on_exit, m_exit)
        
        menu_help = wx.Menu()
        m_about = menu_help.Append(-1, "&About\tF1", "About the demo")
        self.Bind(wx.EVT_MENU, self.on_about, m_about)
        
        self.menubar.Append(menu_file, "&File")
        self.menubar.Append(menu_help, "&Help")
        self.SetMenuBar(self.menubar)
    
    # 2 - Main panel contains the plot, slider bars and output file name entry

	def create_main_panel(self):
    
		self.panel = wx.Panel(self)
		self.dpi = 100
		self.fig = Figure((5.0, 5, 5), dpi = self.dpi)
		self.canvas = FigCanvas(self.panel, -1, self.fig)
		
		self.axes = self.add.add_subplot(111)
    	
    	# 3 - Insert code for tooltip functionality
    	
		self.outputNameBox = wx.TextCtrl(
			self.panel,
			size = (200, -1),
			style = wx.TE_PROCESS_ENTER)
		self.Bind(wx.EVT_TEXT_ENTER, self.on_text_enter, self.outputNameBox)
		
		self.drawbutton = wx.Button(self.panel, -1, "Re-plot")
		self.Bind(wx.EVT_BUTTON, self.on_draw_button, self.drawbutton)
		
		self.save_csv = wx.CheckBox(self.panel, -1,
			"Save CSV",
			style = wx.ALIGN_RIGHT)
		self.Bind(wx.EVT_CHECKBOX, self.on_save_csv, self.save_csv)
		
		self.first_slider_label = wx.StaticText(self.panel, -1, 
			"Scoring cutoff: ")
		self.first_slider_width = wx.Slider(self.panel, -1, 
			value=20, 
			minValue=1,
			maxValue=100,
			style=wx.SL_AUTOTICKS | wx.SL_LABELS)
		self.first_slider_width.SetTickFreq(10, 1)
		self.Bind(wx.EVT_COMMAND_SCROLL_THUMBTRACK, self.on_first_slider_width, self.first_slider_width)
		
		self.toolbar = NavigationToolbar(self.canvas)
		
		self.vbox = wx.BoxSizer(wx.VERTICAL)
		self.vbox.Add(self.canvas, 1, wx.LEFT | wx.TOP | wx.GROW)
		self.vbox.Add(self.toolbar, 0, wx.EXPAND)
		
		self.hbox = wx.BoxSizer(wx.HORIZONTAL)
		flags = wx.ALIGN_LEFT | wx.ALL | wx.ALIGN_CENTER_VERTICAL
		self.hbox.Add(self.outputNameBox, 0, border=3, flag=flags)
		self.hbox.Add(self.drawbutton, 0, border=3, flag=flags)
		self.hbox.Add(self.save_csv, 0, border = 3, flag=flags)
		self.hbox.AddSpacer(10)
		self.hbox.Add(self.first_slider_label, 0, flag=flags)
		self.hbox.Add(self.first_slider_width, 0, border=3, flag=flags)
		
		self.vbox.Add(self.hbox, 0, flag = wx.ALIGN_LEFT | wx.BOTTOM)
		
		self.panel.SetSizer(self.vbox)
		self.vbox.Fit(self)
	
	def create_status_bar(self):
		self.statusbar = self.CreateStatusBar()
	
	def draw_figure(self):
		
		self.axes.clear()
		
		self.axes.set_xlabel('Peptide', fontsize = 20)
		self.axes.set_ylabel('Fraction label', fontsize = 20)
		self.canvas.draw()
    
	def on_save_csv(self, event):
		self.save_csv_option = abs(1-self.save_csv_option)
		print self.save_csv_option
	
	def on_slider_width(self, event):
		self.draw_figure()
	
	def on_draw_button(self, event):
		self.draw_figure()
	
	def on_text_enter(self, event):
		self.draw_figure()
	
	def on_save_plot(self, event):
		file_choices = "PNG (*.png)|*.png"
		
		dlg = wx.FileDialog(
			self, 
			message="Save plot as...",
			defaultDir=os.getcwd(),
			defaultFile="plot.png",
			wildcard=file_choices,
			style=wx.SAVE)
		
		if dlg.ShowModal() == wx.ID_OK:
			path = dlg.GetPath()
			self.canvas.print_figure(path, dpi=self.dpi)
			self.flash_status_message("Saved to %s" % path)
		
	def on_exit(self, event):
		self.Destroy()
						
	def on_about(self, event):
		msg = """ A demo using wxPython with matplotlib:
		
		 * Use the matplotlib navigation bar
		 * Add values to the text box and press Enter (or click "Draw!")
		 * Show or hide the grid
		 * Drag the slider to modify the width of the bars
		 * Save the plot to a file using the File menu
		 * Click on a bar to receive an informative message
		"""
		dlg = wx.MessageDialog(self, msg, "About", wx.OK)
		dlg.ShowModal()
		dlg.Destroy()
		
	def flash_status_message(self, msg, flash_len_ms=1500):
		self.statusbar.SetStatusText(msg)
		self.timeroff = wx.Timer(self)
		self.Bind(
			wx.EVT_TIMER, 
			self.on_flash_status_off, 
			self.timeroff)
		self.timeroff.Start(flash_len_ms, oneShot=True)
	
	def on_flash_status_off(self, event):
		self.statusbar.SetStatusText('')
		
if __name__ == '__main__':
    app = wx.App()
    app.frame = MainFrame()
    app.frame.Show()
    app.MainLoop()