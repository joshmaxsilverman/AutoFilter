import os
import pprint
import random
import wx
import csv
import string
import qMS
import sets
import numpy

# 1 - Code to grab the data, setup the relevant dictionaries

datapath = "/Users/joshsilverman/Dropbox/Research/AutoFilter/gogat_isocsv/gogat40_iso.csv"
datapath = "/Users/joshsilverman/Dropbox/Research/AutoFilter/2-21_iso.csv"

picturepath = "/Users/joshsilverman/Dropbox/Research/AutoFilter/"

data = list( csv.reader( open(datapath, 'rU') ) )
header = data[0]

# 1 - Find the indices for the quantities of interest using list comprehensions

ampu_index = [index for index, item in enumerate(header) if item == "AMP_U"][0]
ampl_index = [index for index, item in enumerate(header) if item == "AMP_L"][0]
isoz_charge_index = [index for index, item in enumerate(header) if item == "isoz_charge"][0]
isopep_index = [index for index, item in enumerate(header) if item == "isopep"][0]
protein_index = [index for index, item in enumerate(header) if item == "protein"][0]
N14ppm_index = [index for index, item in enumerate(header) if item == "ppm_n14"][0]
N15ppm_index = [index for index, item in enumerate(header) if item == "ppm_n15"][0]
UID_index = [index for index, item in enumerate(header) if item == "isofile"][0]
N14abundance_index = [index for index, item in enumerate(header) if item == "abundance_n14"][0]
N15abundance_index = [index for index, item in enumerate(header) if item == "abundance_n15"][0]
missed_index = [index for index, item in enumerate(header) if item == "missed"][0]


# 2 - Declare the peptide list and the protein set

peptide_list = []
protein_set  = set()
peptide_set = set()

data_dictionary = {}
ppm_difference_list = []
N14ppm_list = []
N15ppm_list = []
csv_hold_dict = {}
UID_list = []
missed_cleavage_list = []

# 3 - Loop over data set, collect amplitudes, charge state, peptide sequence, protein id into protein_set

for line in data[1:]:
	
	ampu = float(line[ampu_index])
	ampl = float(line[ampl_index])
	isoz_charge = int(line[isoz_charge_index])
	isopep = line[isopep_index]
	protein = line[protein_index]
	N14ppm = float(line[N14ppm_index])
	N15ppm = float(line[N15ppm_index])
	N14abundance = float(line[N14abundance_index])
	N15abundance = float(line[N15abundance_index])
	missed = float(line[missed_index])
	UID = line[UID_index]	
	ampf = ampu
	
	csv_hold_dict[UID] = line
	UID_list.append(UID)
	
	identifier = [isopep, isoz_charge, protein, ampu/(ampu+ampl)]
	
	data_dictionary[UID] = [ampu/(ampu+ampl), protein, [N14ppm, N15ppm], abs(N14ppm - N15ppm), [N14abundance, N15abundance], missed, ampu, ampl, ampf]
	ppm_difference_list.append( abs(N14ppm - N15ppm) )
	N14ppm_list.append(N14ppm)
	N15ppm_list.append(N15ppm)
	missed_cleavage_list.append(missed)
	
	protein_set.add(protein)
	peptide_set.add(isopep)
	peptide_list.append(identifier)

# ppm_difference_MAD = qMS.MAD(ppm_difference_list)
# N14ppm_MAD = qMS.MAD(N14ppm_list)
# N15ppm_MAD = qMS.MAD(N15ppm_list)
# 
# N14ppm_mean = numpy.mean(N14ppm_list)
# N15ppm_mean = numpy.mean(N15ppm_list)

peptide_MAD_dict = {}
function_scores_dict = {}

# 1 - Filtering by ppm high or low on either peak/difference, filtering by abundance high or low
	
for protein in list(protein_set):
	temp_peptides = filter(lambda item: item[2] == protein, peptide_list)
	temp_fraclabs = map(lambda item: item[3], temp_peptides)
	temp_protein = temp_peptides[0][2]

	peptide_MAD_dict.setdefault(temp_protein , qMS.MAD(temp_fraclabs) )	

protein_index_dict = { item: index for (index, item) in enumerate(list(protein_set)) }

def outputStatsFile(protein_set, filePrefix, numerator, denominator, outputList):
        fileSuffix = "-"
        for i in numerator:
                fileSuffix = fileSuffix+str(i)+"_p_"
        fileSuffix = fileSuffix[:-3]+"_o_"
        for i in denominator:
                fileSuffix = fileSuffix+str(i)+"_p_"
        fileName = filePrefix+fileSuffix[:-3]+".stats"
        outfile = open(fileName, 'w')
        outfile.write("Protein flab loc nval vals\n")
        
        dataByProtein = {}
        for p in protein_set:
        		dataByProtein[p] = {"flab" : -1, "loc" : -1, "nvals" : 0, "vals" : []}
        
        for uid in outputList:
				# Change to grab by hash
                parentProtein = data_dictionary[uid][1]
                u = data_dictionary[uid][6]
                l = data_dictionary[uid][7]
                f = data_dictionary[uid][8]
                valNum = 0
                for i in numerator:
                        if i == "ampu":
                                valNum = valNum + u
                        elif i == "ampl":
                                valNum = valNum + l
                        elif i == "ampf":
                                valNum = valNum + f
                        else:
                                "erorr in num"
                valDen = 0
                for i in denominator:
                        if i == "ampu":
                                valDen = valDen + u
                        elif i == "ampl":
                                valDen = valDen + l
                        elif i == "ampf":
                                valDen = valDen + f
                        else:
                                "erorr in denom"
                value = float(valNum)/float(valDen)
                dataByProtein[parentProtein]["vals"].append(value)
                dataByProtein[parentProtein]["nvals"] = dataByProtein[parentProtein]["nvals"] + 1
        ordered = list(protein_set)
        ordered.sort()
        for p in ordered:
                outString = str(p) + " " + str(numpy.mean(dataByProtein[p]["vals"]))[:6] + " " + str(numpy.std(dataByProtein[p]["vals"]))[:6] + " " + str(dataByProtein[p]["nvals"]) + " "
                for v in dataByProtein[p]["vals"]:
                        outString = outString + str(v)[:6] + " "
                outString = outString + "\n"
                outfile.write(outString)
        outfile.close()



# 2 - The recommended way to use wx with mpl is with the WXAgg backend. 

import matplotlib
matplotlib.use('WXAgg')
from matplotlib.figure import Figure
from matplotlib.backends.backend_wxagg import \
    FigureCanvasWxAgg as FigCanvas, \
    NavigationToolbar2WxAgg as NavigationToolbar


class BarsFrame(wx.Frame):
    """ The main frame of the application
    """
    title = 'Command Station'
    
    def __init__(self):
        wx.Frame.__init__(self, None, -1, self.title)
                
        self.create_menu()
        self.create_status_bar()
        self.create_main_panel()
        
        self.firstRangeBypass.SetValue('0 9')
        self.secondRangeBypass.SetValue('28 43')
        self.thirdRangeBypass.SetValue('25 42')
        self.fourthRangeBypass.SetValue('0 5')

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

    def create_main_panel(self):
        """ Creates the main panel with all the controls on it:
             * mpl canvas 
             * mpl navigation toolbar
             * Control panel for interaction
        """
        self.panel = wx.Panel(self)
        
        self.dpi = 100
        self.fig = Figure((5.0, 4.0), dpi=self.dpi)
        self.canvas = FigCanvas(self.panel, -1, self.fig)
        
        self.axes = self.fig.add_subplot(111)
        
        # Bind the 'pick' event for clicking on one of the bars
        self.canvas.mpl_connect('pick_event', self.on_pick)

        self.first_range_button = wx.Button(self.panel, -1, "PPM difference cutoff")
        self.Bind(wx.EVT_BUTTON, self.on_first_range_button, self.first_range_button)
        
        self.second_range_button = wx.Button(self.panel, -1, "14N PPM cutoff")
        self.Bind(wx.EVT_BUTTON, self.on_second_range_button, self.second_range_button)
        
        self.third_range_button = wx.Button(self.panel, -1, "15N PPM cutoff")
        self.Bind(wx.EVT_BUTTON, self.on_third_range_button, self.third_range_button)
        
        self.fourth_range_button = wx.Button(self.panel, -1, "Missed cleavage cutoff")
        self.Bind(wx.EVT_BUTTON, self.on_fourth_range_button, self.fourth_range_button)
        
############
        
        self.firstRangeBypass = wx.TextCtrl(
            self.panel, 
            size=(200,-1),
            style=wx.TE_PROCESS_ENTER)
        self.Bind(wx.EVT_TEXT_ENTER, self.on_first_range_bypass, self.firstRangeBypass)
        
        self.secondRangeBypass = wx.TextCtrl(
            self.panel, 
            size=(200,-1),
            style=wx.TE_PROCESS_ENTER)
        self.Bind(wx.EVT_TEXT_ENTER, self.on_second_range_bypass, self.secondRangeBypass)

        self.thirdRangeBypass = wx.TextCtrl(
            self.panel, 
            size=(200,-1),
            style=wx.TE_PROCESS_ENTER)
        self.Bind(wx.EVT_TEXT_ENTER, self.on_third_range_bypass, self.thirdRangeBypass)

        self.fourthRangeBypass = wx.TextCtrl(
            self.panel, 
            size=(200,-1),
            style=wx.TE_PROCESS_ENTER)
        self.Bind(wx.EVT_TEXT_ENTER, self.on_fourth_range_bypass, self.fourthRangeBypass)
        
        
        self.exportButton = wx.Button(self.panel, -1, "Export")
        self.Bind(wx.EVT_BUTTON, self.on_export_button, self.exportButton)
        
        self.drawbutton = wx.Button(self.panel, -1, "Draw!")
        self.Bind(wx.EVT_BUTTON, self.on_draw_button, self.drawbutton)

        self.cb_grid = wx.CheckBox(self.panel, -1, 
            "Show Grid",
            style=wx.ALIGN_RIGHT)
        self.Bind(wx.EVT_CHECKBOX, self.on_cb_grid, self.cb_grid)

       # Creates the navigation toolbar
       
        self.toolbar = NavigationToolbar(self.canvas)
        
        # Lays out the various controls
        
        self.vbox = wx.BoxSizer(wx.VERTICAL)
        self.vbox.Add(self.canvas, 1, wx.LEFT | wx.TOP | wx.GROW)
        self.vbox.Add(self.toolbar, 0, wx.EXPAND)

        self.hbox = wx.BoxSizer(wx.HORIZONTAL)
        flags = wx.ALIGN_LEFT | wx.ALL | wx.ALIGN_CENTER_VERTICAL

        self.hbox.Add(self.drawbutton, 0, border=3, flag=flags)
        self.hbox.Add(self.cb_grid, 0, border=3, flag=flags)
        self.hbox.Add(self.exportButton, 0, border = 3, flag = flags) 
      
        self.vbox.Add(self.hbox, 0, flag = wx.ALIGN_LEFT | wx.TOP)
        
        self.vbox.Add((0,0))
        
        # Sliders for setting the various cutoffs
        
        self.sliderBox1 = wx.BoxSizer(wx.HORIZONTAL)
 
 		# The first slider
 		
        self.sliderBox1.Add(self.first_range_button, 0, border=3, flag=flags)
        self.sliderBox1.Add(self.firstRangeBypass, 0, border = 3, flag = flags)

        # The second slider
        self.sliderBox1.Add(self.second_range_button, 0, border=3, flag=flags)
        self.sliderBox1.Add(self.secondRangeBypass, 0, border = 3, flag = flags)
        
        self.vbox.Add(self.sliderBox1, 0, flag = wx.ALIGN_LEFT | wx.TOP)
        
        self.sliderBox2 = wx.BoxSizer(wx.HORIZONTAL)

 		# The third slider
        self.sliderBox2.Add(self.third_range_button, 0, border=3, flag=flags)
        self.sliderBox2.Add(self.thirdRangeBypass, 0, border = 3, flag = flags)

        
 		# The fourth slider
        self.sliderBox2.Add(self.fourth_range_button, 0, border=3, flag=flags)
        self.sliderBox2.Add(self.fourthRangeBypass, 0, border = 3, flag = flags)        

        self.vbox.Add(self.sliderBox2, 0, flag = wx.ALIGN_LEFT | wx.TOP)
        self.panel.SetSizer(self.vbox)
        self.vbox.Fit(self)
        
    def create_status_bar(self):
        self.statusbar = self.CreateStatusBar()

    def draw_figure(self):
    	
    	self.UID_output_list = []
    	
    	self.axes.clear()
    	self.axes.grid(self.cb_grid.IsChecked())
    	
    	for item in data_dictionary.keys():
    	
    		(firstRangeLow, firstRangeHigh) = map(float, self.firstRangeBypass.GetValue().split(' '))
    		(secondRangeLow, secondRangeHigh) = map(float, self.secondRangeBypass.GetValue().split(' '))
    		(thirdRangeLow, thirdRangeHigh) = map(float, self.thirdRangeBypass.GetValue().split(' '))
    		(fourthRangeLow, fourthRangeHigh) = map(float, self.fourthRangeBypass.GetValue().split(' '))
    		
    		# Test the filter conditions
    		con1 = firstRangeLow <= data_dictionary[item][3] <= firstRangeHigh
    		con2 = secondRangeLow <= data_dictionary[item][2][0] <= secondRangeHigh
    		con3 = thirdRangeLow <= data_dictionary[item][2][1] <= thirdRangeHigh
    		con4 = fourthRangeLow <= data_dictionary[item][5] <= fourthRangeHigh
    		
    		# Keep and plot in black peptides that pass the filters, plot in gray those that do not
    		if con1 and con2 and con3 and con4:
    			self.axes.plot( protein_index_dict[data_dictionary[item][1]], data_dictionary[item][0], 'kx')
    			self.UID_output_list.append(item)
    		else:
    			self.axes.plot( protein_index_dict[data_dictionary[item][1]], data_dictionary[item][0], 'rx')
    	self.axes.set_xlim([0, len(protein_index_dict.items())])
    	self.canvas.draw()    
        
    def on_export_button(self, event):
    
    	# Make the list of rejected peptides
    	UID_withheld_list = list(set(UID_list)-set(self.UID_output_list))
    	# Clear the output folder and remake it
    	os.system('rm -r output_withheld_pics')
    	os.system('mkdir output_withheld_pics')
    	# Copy the fit picture of all rejected peptides to the output folder
    	for item in UID_withheld_list:
    		os.system('cp '+picturepath+item+'.fit.png '+'output_withheld_pics')
    		    	
    	# Clear the output folder and remake it
    	os.system('rm -r output_pics')
    	os.system('mkdir output_pics')
    	# Copy the fit picture of all rejected peptides to the output folder
    	for item in self.UID_output_list:
    		os.system('cp '+picturepath+item+'.fit.png '+'output_pics')
    	
    	print 'export complete'

    	# CSV output
    	output_csv = datapath.split('.')[0]+'_filtered.csv'
#     	print output_csv
    	output_csv_open = open(output_csv, 'w')
    	output_csv_open.write(','.join(header))
    	for key in self.UID_output_list:
    		output_csv_open.write(','.join(csv_hold_dict[key])+'\n')
    	output_csv_open.close()
    	
    	outputStatsFile(protein_set, "test", ["ampu"], ["ampu", "ampl"], self.UID_output_list)
    	print self.UID_output_list[0]

    	    
    def on_cb_grid(self, event):
        self.draw_figure()
        
    def on_first_range_bypass(self, event):
    	pass

    def on_second_range_bypass(self, event):
    	pass

    def on_third_range_bypass(self, event):
    	pass

    def on_fourth_range_bypass(self, event):
    	pass
        
    def on_draw_button(self, event):
        self.draw_figure()
    
    def on_first_range_button(self, event):
    	app1 = wx.PySimpleApp()
    	app1.frame = MyFrame("PPM difference cutoff", ppm_difference_list)
    	app1.frame.Show()
    	app.MainLoop() 

    def on_second_range_button(self, event):
    	app1 = wx.PySimpleApp()
    	app1.frame = MyFrame("N14 ppm cutoff", N14ppm_list)
    	app1.frame.Show()
    	app1.MainLoop() 
    	
    def on_third_range_button(self, event):
    	app1 = wx.PySimpleApp()
    	app1.frame = MyFrame("N15 ppm cutoff", N15ppm_list)
    	app1.frame.Show()
    	app1.MainLoop() 
    	
    def on_fourth_range_button(self, event):
    	app1 = wx.PySimpleApp()
    	app1.frame = MyFrame("Missed cleavage cutoff", missed_cleavage_list)
    	app1.frame.Show()
    	app1.MainLoop() 
		
    def on_pick(self, event):
        # The event received here is of the type
        # matplotlib.backend_bases.PickEvent
        #
        # It carries lots of information, of which we're using
        # only a small amount here.
        # 
        box_points = event.artist.get_bbox().get_points()
        msg = "You've clicked on a bar with coords:\n %s" % box_points
        
        dlg = wx.MessageDialog(
            self, 
            msg, 
            "Click!",
            wx.OK | wx.ICON_INFORMATION)

        dlg.ShowModal() 
        dlg.Destroy()        

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
        
# Second window popup for graphs

class MyFrame(wx.Frame):
	title = 'Little douche window'
	
	def __init__(self, passed_title, passed_data):
		self.title = passed_title
		self.data = passed_data
		
		wx.Frame.__init__(self, None, -1, self.title)
		
		self.create_main_panel()
		self.draw_figure(self.data)
	
	def set_title(self, new_title):
		self.title = new_title
		
	def create_main_panel(self):
		""" Creates the main panel with all the controls on it:
		 * mpl canvas 
		 * mpl navigation toolbar
		 * Control panel for interaction
		"""
		self.panel = wx.Panel(self)
		self.dpi = 100
		self.fig = Figure((5.0, 4.0), dpi=self.dpi)
		self.canvas = FigCanvas(self.panel, -1, self.fig)
		
		self.close_button = wx.Button(self.panel, -1, "Close")
		self.Bind(wx.EVT_BUTTON, self.on_close_button, self.close_button)
		
		self.vbox = wx.BoxSizer(wx.VERTICAL)
		self.vbox.Add(self.canvas, 1, wx.LEFT | wx.TOP | wx.GROW)
		self.vbox.Add((0,0))
		
		self.hbox = wx.BoxSizer(wx.HORIZONTAL)
		flags = wx.ALIGN_LEFT | wx.ALL | wx.ALIGN_CENTER_VERTICAL
		self.hbox.Add(self.close_button, 0, border=3, flag=flags)

		self.vbox.Add(self.hbox, 0, flag = wx.ALIGN_CENTER | wx.TOP)
		self.panel.SetSizer(self.vbox)

		self.vbox.Fit(self)
		self.axes = self.fig.add_subplot(111)

		
	def on_close_button(self, event):
		self.Destroy()
		
        
	def draw_figure(self, thing):
		bin_num = min(80, len(list(set(thing))))
		self.axes.hist(thing, bins = bin_num)
		self.canvas.draw()

if __name__ == '__main__':
    app = wx.App()
    app.frame = BarsFrame()
    app.frame.Show()
    app.MainLoop()