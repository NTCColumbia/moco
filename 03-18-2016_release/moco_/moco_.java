
import ij.*;
import ij.measure.ResultsTable;
import ij.process.*;
import ij.gui.*;
import ij.plugin.*;
import ij.plugin.filter.Analyzer;
import ij.io.OpenDialog;


import java.io.IOException;
import java.util.ArrayList;
import java.util.List;

import org.jtransforms.fft.DoubleFFT_2D;


public class moco_ implements PlugIn {
	

	public void run(String arg)  {	
		// Throw error if less than 2 windows are open.
		if (WindowManager.getWindowCount() < 2)
		{
			IJ.error("Open a window containing the template image.");
			return;
		}
		// Obtaining image names
		String names[] = new String[WindowManager.getWindowCount()];
		String templateName;
		String stackName;
		int ids[] = WindowManager.getIDList();
		int sID = 0;
		int tID = 0;
		for (int i = 0; i < ids.length; i++)
		{
			names[i] = WindowManager.getImage(ids[i]).getTitle();
			if (WindowManager.getImage(ids[i]).getImageStackSize() > 1)
			{
				sID = i;
			}
			else
			{
				tID = i;
			}
		}
		// Initialize height, width, ds (downsampling value), and w (maximum width of translation)
		int height = WindowManager.getCurrentImage().getHeight();
		int width  = WindowManager.getCurrentImage().getWidth();
		int ds = 0;
		if (Math.max(height,width) > 256)	{ ds = 1; }
		int w = (int) Math.round((int) Math.min(width,height)/5.0);
		// Choose to generate log file, use existing log file, or none.
		String[] choices = {"Generate log file", "Choose log file", "None"};
		String choice;
		String[] plotChoices = {"Plot RMS", "No plot"};
		String plotChoice;
		// Create Generic Dialog
		GenericDialog gd = new GenericDialog("moco");
		gd.addNumericField("Value of w to use: ", w, 0);
		gd.addNumericField("Downsample_value: ", ds, 0);
		gd.addChoice("Template image",names,names[tID]);
		gd.addChoice("Stack",names,names[sID]);
		gd.addChoice("Log file", choices, choices[2]);
		gd.addChoice("Plot", plotChoices, plotChoices[1]);
		gd.showDialog();
      	if (gd.wasCanceled()) return;
		// Update w and downsample values and template and stack names
		w = (int) gd.getNextNumber();
        ds = (int) gd.getNextNumber();
		templateName = gd.getNextChoice();
		stackName = gd.getNextChoice();
		choice = gd.getNextChoice();
		plotChoice = gd.getNextChoice();
		
		// Throw error if template image is a stack
        if (WindowManager.getImage(templateName).getImageStackSize() > 1)
        {
            IJ.error("Template image must be a single image.");
            return;
        }
		
		// Initialization
        ImagePlus sPlus = WindowManager.getImage(stackName);
        ImagePlus tPlus = WindowManager.getImage(templateName);
        ImageStack stack = sPlus.getImageStack();
        int type = sPlus.getType();
        
        
        // 16-bit image
        if (type == 1)
        {
        
        ShortProcessor proc = (ShortProcessor) tPlus.getProcessor();
        
        int multiplier = 1;
        int size = stack.getSize();
        int rows = proc.getHeight();
        int cols = proc.getWidth();
        int arrSize = ds + 1;
        
        // Template processing
        ShortProcessor[] templates = new ShortProcessor[arrSize];
        templates[0] = (ShortProcessor) proc;
        for (int i = 1; i < arrSize; i++)
        {
            ShortProcessor template = downsample(templates[multiplier-1],rows,cols);
            rows = template.getHeight();
            cols = template.getWidth();
            multiplier++;
            templates[multiplier-1] = template;
        }
        double[][] tVals = computeTVals(templates[arrSize-1],rows,cols);
        List<double[][]> tNew = computeT(tVals,cols,rows);
        w = (int) Math.floor(Math.min(Math.min(w+1, rows), cols)/(Math.pow(2, multiplier-1))) - 1;
        if (w < 0)
        {
            w = 1;
        }
        double[][] tFFT = templateFFT(tVals,rows,cols,w);
        DoubleFFT_2D fft2D = new DoubleFFT_2D(cols+w,rows+w);
        
        boolean ce = false;
        boolean decision = false;
        ResultsTable rt = new ResultsTable();
        ImageStack newStack = new ImageStack(width,height);
		
		
		// Generate log file choice
		if (choice.equals(choices[0]))
		{
			Analyzer.setResultsTable(rt);
			
			for (int t=1; t <= size; t++)
	        {
	        	ShortProcessor image = (ShortProcessor) stack.getProcessor(t);
				ShortProcessor[] imageStack = downsample(image,height,width,arrSize);
	        	int[] xy = templateTranslation(imageStack[arrSize-1],imageStack[arrSize-1].getHeight(),
	        		imageStack[arrSize-1].getWidth(),tFFT,w,fft2D, tNew.get(0),tNew.get(1),tNew.get(2),tNew.get(3));
	        	
	        	for (int i=arrSize-1-1; i>=0; i--)
	        	{
	        		xy[0] = 2*xy[0];
	        		xy[1] = 2*xy[1];
	        		xy = moveByOne2(xy,imageStack[i],templates[i]);
	        	}
	        	
				ShortProcessor newImage = applyAffine(xy,image,height,width);  
				newStack.addSlice(newImage);
				System.out.println("Index " + t + ": " + xy[0] + ", " + xy[1]);
				rt.incrementCounter();
				rt.addValue("x", xy[0]);
				rt.addValue("y", xy[1]);
				IJ.showStatus("Registering image " + t + "/" + size);
				
				// Throw a warning if the translation terms are large.
				if ((Math.abs(xy[0]) >= 30 || Math.abs(xy[1]) >= 30) && decision == false)
                {
                    GenericDialog gdce = new GenericDialog("Contrast Enhancement");
                    gdce.addMessage("The translation term is very large. Contrast enhancement may ensure that the results are accurate. \n Would you like to apply contrast enhancement?");
                    gdce.enableYesNoCancel();
                    gdce.hideCancelButton();
                    gdce.showDialog();
                    if (gdce.wasOKed())
                    {
                        ce = true;
                        break;
                    }
                    else
                    {
                        decision = true;
                    }
                }
				
	        }
			
		}
		// Choose log file
		else if (choice.equals(choices[1]))
		{
			OpenDialog dialog = new OpenDialog("Choose the file to open.");
			String path = dialog.getPath();
			try {
				rt = ResultsTable.open(path);
			} catch (IOException e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
			}
			rt.show("Results");
			for (int t=1; t <= size; t++)
	        {
	        	ShortProcessor image = (ShortProcessor) stack.getProcessor(t);
	        	int[] xy = new int[2];
	        	xy[0] = (int) rt.getValueAsDouble(0, t-1);
	        	xy[1] = (int) rt.getValueAsDouble(1, t-1);
				ShortProcessor newImage = applyAffine(xy,image,height,width);  
				newStack.addSlice(newImage);
				System.out.println("Index " + t + ": " + xy[0] + ", " + xy[1]);
				IJ.showStatus("Registering image " + t + "/" + size);
	        }
		}
		// Neither choose nor generate log file
		else if (choice.equals(choices[2]))
		{
			for (int t=1; t <= size; t++)
	        {
	        	ShortProcessor image = (ShortProcessor) stack.getProcessor(t);
				ShortProcessor[] imageStack = downsample(image,height,width,arrSize);
	        	int[] xy = templateTranslation(imageStack[arrSize-1],imageStack[arrSize-1].getHeight(),
	        		imageStack[arrSize-1].getWidth(),tFFT,w,fft2D, tNew.get(0),tNew.get(1),tNew.get(2),tNew.get(3));
	        	
	        	for (int i=arrSize-1-1; i>=0; i--)
	        	{
	        		xy[0] = 2*xy[0];
	        		xy[1] = 2*xy[1];
	        		xy = moveByOne2(xy,imageStack[i],templates[i]);
	        	}
	        	
				ShortProcessor newImage = applyAffine(xy,image,height,width);  
				newStack.addSlice(newImage);
				System.out.println("Index " + t + ": " + xy[0] + ", " + xy[1]);
				IJ.showStatus("Registering image " + t + "/" + size);
				if ((Math.abs(xy[0]) >= 30 || Math.abs(xy[1]) >= 30) && decision == false)
                {
                    GenericDialog gdce = new GenericDialog("Contrast Enhancement");
                    gdce.addMessage("The translation term is very large. Contrast enhancement may ensure that the results are accurate. \n Would you like to apply contrast enhancement?");
                    gdce.enableYesNoCancel();
                    gdce.hideCancelButton();
                    gdce.showDialog();
                    if (gdce.wasOKed())
                    {
                        ce = true;
                        break;
                    }
                    else
                    {
                        decision = true;
                    }
                }
	        }
		}
        
		
		
		// Apply contrast enhancement before performing the moco algorithm
        if (ce == true)
        {
            IJ.selectWindow(stackName);
            IJ.run("Enhance Contrast...", "saturated=0.3 normalize process_all");
            IJ.selectWindow(templateName);
            IJ.run("Enhance Contrast...", "saturated=0.3 normalize");
            // Template processing (again)
            sPlus = WindowManager.getImage(stackName);
            tPlus = WindowManager.getImage(templateName);
            stack = sPlus.getImageStack();
            proc = (ShortProcessor) tPlus.getProcessor();
            templates = new ShortProcessor[arrSize];
            templates[0] = (ShortProcessor) proc;
            multiplier = 1;
            size = stack.getSize();
            type = sPlus.getType();
            rows = proc.getHeight();
            cols = proc.getWidth();
            arrSize = ds + 1;
            for (int i = 1; i < arrSize; i++)
            {
                ShortProcessor template = downsample(templates[multiplier-1],rows,cols);
                rows = template.getHeight();
                cols = template.getWidth();
                multiplier++;
                templates[multiplier-1] = template;
            }
            tVals = computeTVals(templates[arrSize-1],rows,cols);
            tNew = computeT(tVals,cols,rows);

            tFFT = templateFFT(tVals,rows,cols,w);
            fft2D = new DoubleFFT_2D(cols+w,rows+w);
            
            newStack = new ImageStack(width,height);
            if (choice.equals(choices[0]))
            {
                rt = new ResultsTable();
                Analyzer.setResultsTable(rt);
                for (int t=1; t <= size; t++)
                {
                    ShortProcessor image = (ShortProcessor) stack.getProcessor(t);
                    ShortProcessor[] imageStack = downsample(image,height,width,arrSize);
                    int[] xy = templateTranslation(imageStack[arrSize-1],imageStack[arrSize-1].getHeight(),
                                                   imageStack[arrSize-1].getWidth(),tFFT,w,fft2D, tNew.get(0),tNew.get(1),tNew.get(2),tNew.get(3));
                    
                    for (int i=arrSize-1-1; i>=0; i--)
                    {
                        xy[0] = 2*xy[0];
                        xy[1] = 2*xy[1];
                        xy = moveByOne2(xy,imageStack[i],templates[i]);
                    }
                    
                    ShortProcessor newImage = applyAffine(xy,image,height,width);
                    newStack.addSlice(newImage);
                    System.out.println("Index " + t + ": " + xy[0] + ", " + xy[1]);
                    rt.incrementCounter();
                    rt.addValue("x", xy[0]);
                    rt.addValue("y", xy[1]);
                    IJ.showStatus("Registering image " + t + "/" + size);
                }
            }
            else if (choice.equals(choices[2]))
            {
                for (int t=1; t <= size; t++)
                {
                    ShortProcessor image = (ShortProcessor) stack.getProcessor(t);
                    ShortProcessor[] imageStack = downsample(image,height,width,arrSize);
                    int[] xy = templateTranslation(imageStack[arrSize-1],imageStack[arrSize-1].getHeight(),
                                                   imageStack[arrSize-1].getWidth(),tFFT,w,fft2D, tNew.get(0),tNew.get(1),tNew.get(2),tNew.get(3));
                    
                    for (int i=arrSize-1-1; i>=0; i--)
                    {
                        xy[0] = 2*xy[0];
                        xy[1] = 2*xy[1];
                        xy = moveByOne2(xy,imageStack[i],templates[i]);
                    }
                    
                    ShortProcessor newImage = applyAffine(xy,image,height,width);
                    newStack.addSlice(newImage);
                    System.out.println("Index " + t + ": " + xy[0] + ", " + xy[1]);
                    IJ.showStatus("Registering image " + t + "/" + size);
                }
            }
        }
		
		// Finalize
        ImagePlus newImages = new ImagePlus("New Stack", newStack);
        StackWindow window = new StackWindow(newImages);
        if (choice.equals(choices[0]))
        {
            rt.show("Results");
        }
		if (plotChoice.equals(plotChoices[0]) && choice.equals(choices[0]))
		{
			double[] rms = new double[size];
			double[] frm = new double[size];
			for (int t=0; t < size; t++)
			{
				rms[t] = Math.sqrt(rt.getValueAsDouble(0,t)*rt.getValueAsDouble(0,t) + rt.getValueAsDouble(1,t)*rt.getValueAsDouble(1,t));
				frm[t] = t;
			}
			Plot rmsPlot = new Plot("Errors","Frame","RMS",frm,rms);
			rmsPlot.show();
		}
        System.out.println("Image stablizations are done!");
        }
        
        // 8-bit image
        else if (type == 0)
        {
            ByteProcessor proc = (ByteProcessor) tPlus.getProcessor();
            
            int multiplier = 1;
            int size = stack.getSize();
            int rows = proc.getHeight();
            int cols = proc.getWidth();
            int arrSize = ds + 1;
            
            // Template processing
            ByteProcessor[] templates = new ByteProcessor[arrSize];
            templates[0] = (ByteProcessor) proc;
            for (int i = 1; i < arrSize; i++)
            {
                ByteProcessor template = downsample(templates[multiplier-1],rows,cols);
                rows = template.getHeight();
                cols = template.getWidth();
                multiplier++;
                templates[multiplier-1] = template;
            }
            double[][] tVals = computeTVals(templates[arrSize-1],rows,cols);
            List<double[][]> tNew = computeT(tVals,cols,rows);
            w = (int) Math.floor(Math.min(Math.min(w+1, rows), cols)/(Math.pow(2, multiplier-1))) - 1;
            if (w < 0)
            {
                w = 1;
            }
            double[][] tFFT = templateFFT(tVals,rows,cols,w);
            DoubleFFT_2D fft2D = new DoubleFFT_2D(cols+w,rows+w);
            
            boolean ce = false;
            boolean decision = false;
            ResultsTable rt = new ResultsTable();
            ImageStack newStack = new ImageStack(width,height);
    		
    		
    		// Generate log file choice
    		if (choice.equals(choices[0]))
    		{
    			Analyzer.setResultsTable(rt);
    			
    			for (int t=1; t <= size; t++)
    	        {
    	        	ByteProcessor image = (ByteProcessor) stack.getProcessor(t);
    				ByteProcessor[] imageStack = downsample(image,height,width,arrSize);
    	        	int[] xy = templateTranslation(imageStack[arrSize-1],imageStack[arrSize-1].getHeight(),
    	        		imageStack[arrSize-1].getWidth(),tFFT,w,fft2D, tNew.get(0),tNew.get(1),tNew.get(2),tNew.get(3));
    	        	
    	        	for (int i=arrSize-1-1; i>=0; i--)
    	        	{
    	        		xy[0] = 2*xy[0];
    	        		xy[1] = 2*xy[1];
    	        		xy = moveByOne2(xy,imageStack[i],templates[i]);
    	        	}
    	        	
    				ByteProcessor newImage = applyAffine(xy,image,height,width);  
    				newStack.addSlice(newImage);
    				System.out.println("Index " + t + ": " + xy[0] + ", " + xy[1]);
    				rt.incrementCounter();
    				rt.addValue("x", xy[0]);
    				rt.addValue("y", xy[1]);
    				IJ.showStatus("Registering image " + t + "/" + size);
    				
    				// Throw a warning if the translation terms are large.
    				if ((Math.abs(xy[0]) >= 30 || Math.abs(xy[1]) >= 30) && decision == false)
                    {
                        GenericDialog gdce = new GenericDialog("Contrast Enhancement");
                        gdce.addMessage("The translation term is very large. Contrast enhancement may ensure that the results are accurate. \n Would you like to apply contrast enhancement?");
                        gdce.enableYesNoCancel();
                        gdce.hideCancelButton();
                        gdce.showDialog();
                        if (gdce.wasOKed())
                        {
                            ce = true;
                            break;
                        }
                        else
                        {
                            decision = true;
                        }
                    }
    				
    	        }
    			
    		}
    		// Choose log file
    		else if (choice.equals(choices[1]))
    		{
    			OpenDialog dialog = new OpenDialog("Choose the file to open.");
    			String path = dialog.getPath();
    			try {
    				rt = ResultsTable.open(path);
    			} catch (IOException e) {
    				// TODO Auto-generated catch block
    				e.printStackTrace();
    			}
    			rt.show("Results");
    			for (int t=1; t <= size; t++)
    	        {
    	        	ByteProcessor image = (ByteProcessor) stack.getProcessor(t);
    	        	int[] xy = new int[2];
    	        	xy[0] = (int) rt.getValueAsDouble(0, t-1);
    	        	xy[1] = (int) rt.getValueAsDouble(1, t-1);
    				ByteProcessor newImage = applyAffine(xy,image,height,width);  
    				newStack.addSlice(newImage);
    				System.out.println("Index " + t + ": " + xy[0] + ", " + xy[1]);
    				IJ.showStatus("Registering image " + t + "/" + size);
    	        }
    		}
    		// Neither choose nor generate log file
    		else if (choice.equals(choices[2]))
    		{
    			for (int t=1; t <= size; t++)
    	        {
    	        	ByteProcessor image = (ByteProcessor) stack.getProcessor(t);
    				ByteProcessor[] imageStack = downsample(image,height,width,arrSize);
    	        	int[] xy = templateTranslation(imageStack[arrSize-1],imageStack[arrSize-1].getHeight(),
    	        		imageStack[arrSize-1].getWidth(),tFFT,w,fft2D, tNew.get(0),tNew.get(1),tNew.get(2),tNew.get(3));
    	        	
    	        	for (int i=arrSize-1-1; i>=0; i--)
    	        	{
    	        		xy[0] = 2*xy[0];
    	        		xy[1] = 2*xy[1];
    	        		xy = moveByOne2(xy,imageStack[i],templates[i]);
    	        	}
    	        	
    				ByteProcessor newImage = applyAffine(xy,image,height,width);  
    				newStack.addSlice(newImage);
    				System.out.println("Index " + t + ": " + xy[0] + ", " + xy[1]);
    				IJ.showStatus("Registering image " + t + "/" + size);
    				if ((Math.abs(xy[0]) >= 30 || Math.abs(xy[1]) >= 30) && decision == false)
                    {
                        GenericDialog gdce = new GenericDialog("Contrast Enhancement");
                        gdce.addMessage("The translation term is very large. Contrast enhancement may ensure that the results are accurate. \n Would you like to apply contrast enhancement?");
                        gdce.enableYesNoCancel();
                        gdce.hideCancelButton();
                        gdce.showDialog();
                        if (gdce.wasOKed())
                        {
                            ce = true;
                            break;
                        }
                        else
                        {
                            decision = true;
                        }
                    }
    	        }
    		}
            
    		
    		
    		// Apply contrast enhancement before performing the moco algorithm
            if (ce == true)
            {
                IJ.selectWindow(stackName);
                IJ.run("Enhance Contrast...", "saturated=0.3 normalize process_all");
                IJ.selectWindow(templateName);
                IJ.run("Enhance Contrast...", "saturated=0.3 normalize");
                // Template processing (again)
                sPlus = WindowManager.getImage(stackName);
                tPlus = WindowManager.getImage(templateName);
                stack = sPlus.getImageStack();
                proc = (ByteProcessor) tPlus.getProcessor();
                templates = new ByteProcessor[arrSize];
                templates[0] = (ByteProcessor) proc;
                multiplier = 1;
                size = stack.getSize();
                type = sPlus.getType();
                rows = proc.getHeight();
                cols = proc.getWidth();
                arrSize = ds + 1;
                for (int i = 1; i < arrSize; i++)
                {
                    ByteProcessor template = downsample(templates[multiplier-1],rows,cols);
                    rows = template.getHeight();
                    cols = template.getWidth();
                    multiplier++;
                    templates[multiplier-1] = template;
                }
                tVals = computeTVals(templates[arrSize-1],rows,cols);
                tNew = computeT(tVals,cols,rows);

                tFFT = templateFFT(tVals,rows,cols,w);
                fft2D = new DoubleFFT_2D(cols+w,rows+w);
                
                newStack = new ImageStack(width,height);
                if (choice.equals(choices[0]))
                {
                    rt = new ResultsTable();
                    Analyzer.setResultsTable(rt);
                    for (int t=1; t <= size; t++)
                    {
                        ByteProcessor image = (ByteProcessor) stack.getProcessor(t);
                        ByteProcessor[] imageStack = downsample(image,height,width,arrSize);
                        int[] xy = templateTranslation(imageStack[arrSize-1],imageStack[arrSize-1].getHeight(),
                                                       imageStack[arrSize-1].getWidth(),tFFT,w,fft2D, tNew.get(0),tNew.get(1),tNew.get(2),tNew.get(3));
                        
                        for (int i=arrSize-1-1; i>=0; i--)
                        {
                            xy[0] = 2*xy[0];
                            xy[1] = 2*xy[1];
                            xy = moveByOne2(xy,imageStack[i],templates[i]);
                        }
                        
                        ByteProcessor newImage = applyAffine(xy,image,height,width);
                        newStack.addSlice(newImage);
                        System.out.println("Index " + t + ": " + xy[0] + ", " + xy[1]);
                        rt.incrementCounter();
                        rt.addValue("x", xy[0]);
                        rt.addValue("y", xy[1]);
                        IJ.showStatus("Registering image " + t + "/" + size);
                    }
                }
                else if (choice.equals(choices[2]))
                {
                    for (int t=1; t <= size; t++)
                    {
                        ByteProcessor image = (ByteProcessor) stack.getProcessor(t);
                        ByteProcessor[] imageStack = downsample(image,height,width,arrSize);
                        int[] xy = templateTranslation(imageStack[arrSize-1],imageStack[arrSize-1].getHeight(),
                                                       imageStack[arrSize-1].getWidth(),tFFT,w,fft2D, tNew.get(0),tNew.get(1),tNew.get(2),tNew.get(3));
                        
                        for (int i=arrSize-1-1; i>=0; i--)
                        {
                            xy[0] = 2*xy[0];
                            xy[1] = 2*xy[1];
                            xy = moveByOne2(xy,imageStack[i],templates[i]);
                        }
                        
                        ByteProcessor newImage = applyAffine(xy,image,height,width);
                        newStack.addSlice(newImage);
                        System.out.println("Index " + t + ": " + xy[0] + ", " + xy[1]);
                        IJ.showStatus("Registering image " + t + "/" + size);
                    }
                }
            }
    		
    		// Finalize
            ImagePlus newImages = new ImagePlus("New Stack", newStack);
            StackWindow window = new StackWindow(newImages);
            if (choice.equals(choices[0]))
            {
                rt.show("Results");
            }
    		if (plotChoice.equals(plotChoices[0]) && choice.equals(choices[0]))
    		{
    			double[] rms = new double[size];
    			double[] frm = new double[size];
    			for (int t=0; t < size; t++)
    			{
    				rms[t] = Math.sqrt(rt.getValueAsDouble(0,t)*rt.getValueAsDouble(0,t) + rt.getValueAsDouble(1,t)*rt.getValueAsDouble(1,t));
    				frm[t] = t;
    			}
    			Plot rmsPlot = new Plot("Errors","Frame","RMS",frm,rms);
    			rmsPlot.show();
    		}
            System.out.println("Image stablizations are done!");

        }
        
        else
        {
        	IJ.error("The stack must be of type 8-bit or 16-bit.");
            return;
        }
		
	}

	public static ShortProcessor downsample(ShortProcessor srcImage, int rows, int cols)
	{
		int dstRows = (int) Math.floor(rows/2);
		int dstCols = (int) Math.floor(cols/2);
		int rowCount = 0;
		int colCount = 0;
		ShortProcessor dstImage = new ShortProcessor(dstCols, dstRows);
		while (rowCount < dstRows)
		{
			colCount = 0;
			while (colCount < dstCols)
			{
				dstImage.set(colCount+dstCols*rowCount,(int) (.25*(srcImage.get(2*colCount+2*rowCount*cols)
					+ srcImage.get(2*colCount+1+2*rowCount*cols) + srcImage.get(2*colCount + (2*rowCount+1)*cols)
					+ srcImage.get(2*colCount+1 + (2*rowCount+1)*cols))));
				
				colCount++;
			}
			rowCount++;
		}
		return dstImage;
	}
	
	public static ByteProcessor downsample(ByteProcessor srcImage, int rows, int cols)
	{
		int dstRows = (int) Math.floor(rows/2);
		int dstCols = (int) Math.floor(cols/2);
		int rowCount = 0;
		int colCount = 0;
		ByteProcessor dstImage = new ByteProcessor(dstCols, dstRows);
		while (rowCount < dstRows)
		{
			colCount = 0;
			while (colCount < dstCols)
			{
				dstImage.set(colCount+dstCols*rowCount,(int) (.25*(srcImage.get(2*colCount+2*rowCount*cols)
					+ srcImage.get(2*colCount+1+2*rowCount*cols) + srcImage.get(2*colCount + (2*rowCount+1)*cols)
					+ srcImage.get(2*colCount+1 + (2*rowCount+1)*cols))));
				
				colCount++;
			}
			rowCount++;
		}
		return dstImage;
	}

	public static double[][] computeTVals(ShortProcessor template, int rows, int cols)
	{
		double[] ms = meanAndStdDev(template,rows,cols);
		double[][] tVals = new double[cols][rows];
		for (int r = 0; r < rows; r++)
		{
			for (int c = 0; c < cols; c++)
			{
				tVals[c][r] = (template.get(c + cols*r) - ms[0])/(ms[1]*Math.sqrt(2));
			}
		}
		return tVals;
		
	}
	
	public static double[][] computeTVals(ByteProcessor template, int rows, int cols)
	{
		double[] ms = meanAndStdDev(template,rows,cols);
		double[][] tVals = new double[cols][rows];
		for (int r = 0; r < rows; r++)
		{
			for (int c = 0; c < cols; c++)
			{
				tVals[c][r] = (template.get(c + cols*r) - ms[0])/(ms[1]*Math.sqrt(2));
			}
		}
		return tVals;
		
	}

	public static List<double[][]> computeT(double[][] tVals, int cols, int rows)
	{
		List<double[][]> T = new ArrayList<double[][]>();

		
		
		double[][] t00Vals = new double[cols][rows];
		double[][] t10Vals = new double[cols][rows];
		double[][] t01Vals = new double[cols][rows];
		double[][] t11Vals = new double[cols][rows];
		t00Vals[0][0] = tVals[0][0]*tVals[0][0];
		t10Vals[0][0] = tVals[0][rows-1]*tVals[0][rows-1];
		t01Vals[0][0] = tVals[cols-1][0]*tVals[cols-1][0];
		t11Vals[0][0] = tVals[cols-1][rows-1]*tVals[cols-1][rows-1];
		for (int i=1; i < rows; i++)
		{	  
			t00Vals[0][i] = t00Vals[0][i-1]+tVals[0][i]*tVals[0][i];
			t10Vals[0][i] = t10Vals[0][i-1]+tVals[0][rows-i-1]*tVals[0][rows-i-1];
			t01Vals[0][i] = t01Vals[0][i-1]+tVals[cols-1][i]*tVals[cols-1][i];
			t11Vals[0][i] = t11Vals[0][i-1]+tVals[cols-1][rows-i-1]*tVals[cols-1][rows-i-1];
		}
		for (int j=1; j< cols; j++)
		{	  
			t00Vals[j][0] = t00Vals[j-1][0]+tVals[j][0]*tVals[j][0];
			t10Vals[j][0] = t10Vals[j-1][0]+tVals[j][rows-1]*tVals[j][rows-1];
			t01Vals[j][0] = t01Vals[j-1][0]+tVals[cols-j-1][0]*tVals[cols-j-1][0];
			t11Vals[j][0] = t11Vals[j-1][0]+tVals[cols-j-1][rows-1]*tVals[cols-j-1][rows-1];
			for (int i = 1; i < rows; i++)
			{
				t00Vals[j][i] = t00Vals[j][i-1]+t00Vals[j-1][i]-t00Vals[j-1][i-1]+tVals[j][i]*tVals[j][i];
				t10Vals[j][i] = t10Vals[j][i-1]+t10Vals[j-1][i]-t10Vals[j-1][i-1]+tVals[j][rows-i-1]*tVals[j][rows-i-1];
				t01Vals[j][i] = t01Vals[j][i-1]+t01Vals[j-1][i]-t01Vals[j-1][i-1]+tVals[cols-j-1][i]*tVals[cols-j-1][i];
				t11Vals[j][i] = t11Vals[j][i-1]+t11Vals[j-1][i]-t11Vals[j-1][i-1]+tVals[cols-j-1][rows-i-1]*tVals[cols-j-1][rows-i-1];
			} 
		}
		
		T.add(t00Vals);
		T.add(t10Vals);
		T.add(t01Vals);
		T.add(t11Vals);
		
		return T;	  	  
	}
	
	
	

	public static ShortProcessor[] downsample(ShortProcessor image, int rows, int cols, int arrSize)
	{
		ShortProcessor[] imageStack = new ShortProcessor[arrSize];
		imageStack[0] = image;
		for (int i=1; i<arrSize; i++)
		{
			imageStack[i] = downsample(imageStack[i-1],rows,cols);
			rows = imageStack[i].getHeight();
			cols = imageStack[i].getWidth();
		}
		
		return imageStack;
	}
	
	public static ByteProcessor[] downsample(ByteProcessor image, int rows, int cols, int arrSize)
	{
		ByteProcessor[] imageStack = new ByteProcessor[arrSize];
		imageStack[0] = image;
		for (int i=1; i<arrSize; i++)
		{
			imageStack[i] = downsample(imageStack[i-1],rows,cols);
			rows = imageStack[i].getHeight();
			cols = imageStack[i].getWidth();
		}
		
		return imageStack;
	}

	
	public static double[][] templateFFT(double[][] oldTVals, int rows, int cols, int w)
	{
		double[][] newTVals = new double[cols][rows];
		for (int r=0; r<rows; r++)
		{
			for (int c=0; c<cols; c++)
			{
				newTVals[c][r] = oldTVals[c][r];
			}
		}
		
		
		  for (int r=0; r<rows; r++)
		  {
			  for (int c=0; c<cols/2; c++)
			  {
				  double x = newTVals[cols-1-c][r];
				  newTVals[cols-1-c][r] = newTVals[c][r];
				  newTVals[c][r] = x;
			  }
		  }
		  // flip each row
		  for (int c=0; c<cols; c++)
		  {
			  for (int r=0; r<rows/2; r++)
			  {
				  double y = newTVals[c][rows-1-r];
				  newTVals[c][rows-1-r] = newTVals[c][r];
				  newTVals[c][r] = y;
			  }
		  }
		  
		  // place values in array for FFT 
		  double[][] tFFT = new double[cols+w][2*(rows+w)];
		  for (int r = 0; r<rows; r++)
		  {
			  for (int c = 0; c< cols; c++)
			  {
				  tFFT[c][r] = newTVals[c][r];
			  }
		  }
		  
		  // Take fft of temp
		  DoubleFFT_2D fft2D = new DoubleFFT_2D(cols+w,rows+w);
		  fft2D.realForwardFull(tFFT);
		  return tFFT;
		  


	}
	

	
	public static int[] templateTranslation(ShortProcessor newImg, int rows, int cols, double[][] tFFT, int w, DoubleFFT_2D fft2D,
			double[][] t00Vals, double[][] t10Vals, double[][] t01Vals, double[][] t11Vals)
	{
		double[] msa = meanAndStdDev(newImg,rows,cols);
		double[][] aVals = new double[cols+w][2*(rows+w)];
		for (int r=0; r<rows; r++)
		{
			for (int c=0; c<cols; c++)
			{
				aVals[c][r] = (newImg.get(c + r*cols) - msa[0])/(msa[1]*Math.sqrt(2));
			}
		}
		
		double[][] a00Vals = new double[cols][rows];
		double[][] a10Vals = new double[cols][rows];
		double[][] a01Vals = new double[cols][rows];
		double[][] a11Vals = new double[cols][rows];
		a00Vals[0][0] = aVals[0][0]*aVals[0][0];
		a10Vals[0][0] = aVals[0][rows-1]*aVals[0][rows-1];
		a01Vals[0][0] = aVals[cols-1][0]*aVals[cols-1][0];
		a11Vals[0][0] = aVals[cols-1][rows-1]*aVals[cols-1][rows-1];
		for (int i=1; i < rows; i++)
		{
			a00Vals[0][i] = a00Vals[0][i-1]+aVals[0][i]*aVals[0][i];
			a10Vals[0][i] = a10Vals[0][i-1]+aVals[0][rows-i-1]*aVals[0][rows-i-1];
			a01Vals[0][i] = a01Vals[0][i-1]+aVals[cols-1][i]*aVals[cols-1][i];
			a11Vals[0][i] = a11Vals[0][i-1]+aVals[cols-1][rows-i-1]*aVals[cols-1][rows-i-1];  
		}
		for (int j=1; j< cols; j++)
   		{	  
			a00Vals[j][0] = a00Vals[j-1][0]+aVals[j][0]*aVals[j][0];
			a10Vals[j][0] = a10Vals[j-1][0]+aVals[j][rows-1]*aVals[j][rows-1];
			a01Vals[j][0] = a01Vals[j-1][0]+aVals[cols-j-1][0]*aVals[cols-j-1][0];
			a11Vals[j][0] = a11Vals[j-1][0]+aVals[cols-j-1][rows-1]*aVals[cols-j-1][rows-1];  
			for (int i = 1; i < rows; i++)
			{
				a00Vals[j][i] = a00Vals[j][i-1]+a00Vals[j-1][i]-a00Vals[j-1][i-1]+aVals[j][i]*aVals[j][i];
				a10Vals[j][i] = a10Vals[j][i-1]+a10Vals[j-1][i]-a10Vals[j-1][i-1]+aVals[j][rows-i-1]*aVals[j][rows-i-1];
				a01Vals[j][i] = a01Vals[j][i-1]+a01Vals[j-1][i]-a01Vals[j-1][i-1]+aVals[cols-j-1][i]*aVals[cols-j-1][i];
				a11Vals[j][i] = a11Vals[j][i-1]+a11Vals[j-1][i]-a11Vals[j-1][i-1]+aVals[cols-j-1][rows-i-1]*aVals[cols-j-1][rows-i-1];  
			} 
		}
		double[][] b = new double[cols+w][rows+w];
		int rowCount = 0; int colCount = 0;

		  

		fft2D.realForwardFull(aVals);
		double[][] out = multFFT(tFFT, aVals, cols+w, 2*(rows+w));	// multiply fft_temp and fft_a
		fft2D.complexInverse(out, true);
		colCount = 0;
		for (int c=0; c<cols+w; c++)
		{
			rowCount = 0;
			for (int r=0; r<2*(rows+w); r+=2)
			{
				b[colCount][rowCount] = out[c][r];
				rowCount++;
			}
			colCount++;
		}
		int[] xy = new int[2];
		double minner = Double.POSITIVE_INFINITY;
		double minner_new = 0;
		for (int r = rows-w-1; r < rows+w; r++)
		{
			for (int c = cols-w-1; c < cols+w; c++)
			{
				if (r < rows-1 && c < cols-1)
				{
					minner_new = (a11Vals[c][r] + t00Vals[c][r] - 2*b[2*(cols-1)-c][2*(rows-1)-r])/(r*c);
				}
				else if (r >= rows-1 && c < cols-1)
				{
				    minner_new = (a01Vals[c][2*(rows-1)-r] + t10Vals[c][2*(rows-1)-r] - 2*b[2*(cols-1)-c][2*(rows-1)-r])/((2*rows-r)*c);
				}
				else if (r < rows-1 && c >= cols-1)
				{
				  	minner_new = (a10Vals[2*(cols-1)-c][r] + t01Vals[2*(cols-1)-c][r] - 2*b[2*(cols-1)-c][2*(rows-1)-r])/(r*(2*cols-c));
				}
				else 
				{
				   	minner_new = (a00Vals[2*(cols-1)-c][2*(rows-1)-r] + t11Vals[2*(cols-1)-c][2*(rows-1)-r] - 2*b[2*(cols-1)-c][2*(rows-1)-r])/((2*rows-r)*(2*cols-c));
				}

				if (minner_new < minner)
				{
					xy[0] = c; xy[1] = r; 
					minner = minner_new;
				}
			}
		}
		xy[0] = xy[0] - (cols-1);
		xy[1] = xy[1] - (rows-1);
		  
		return xy;
	}
	
	public static int[] templateTranslation(ByteProcessor newImg, int rows, int cols, double[][] tFFT, int w, DoubleFFT_2D fft2D,
			double[][] t00Vals, double[][] t10Vals, double[][] t01Vals, double[][] t11Vals)
	{
		double[] msa = meanAndStdDev(newImg,rows,cols);
		double[][] aVals = new double[cols+w][2*(rows+w)];
		for (int r=0; r<rows; r++)
		{
			for (int c=0; c<cols; c++)
			{
				aVals[c][r] = (newImg.get(c + r*cols) - msa[0])/(msa[1]*Math.sqrt(2));
			}
		}
		
		double[][] a00Vals = new double[cols][rows];
		double[][] a10Vals = new double[cols][rows];
		double[][] a01Vals = new double[cols][rows];
		double[][] a11Vals = new double[cols][rows];
		a00Vals[0][0] = aVals[0][0]*aVals[0][0];
		a10Vals[0][0] = aVals[0][rows-1]*aVals[0][rows-1];
		a01Vals[0][0] = aVals[cols-1][0]*aVals[cols-1][0];
		a11Vals[0][0] = aVals[cols-1][rows-1]*aVals[cols-1][rows-1];
		for (int i=1; i < rows; i++)
		{
			a00Vals[0][i] = a00Vals[0][i-1]+aVals[0][i]*aVals[0][i];
			a10Vals[0][i] = a10Vals[0][i-1]+aVals[0][rows-i-1]*aVals[0][rows-i-1];
			a01Vals[0][i] = a01Vals[0][i-1]+aVals[cols-1][i]*aVals[cols-1][i];
			a11Vals[0][i] = a11Vals[0][i-1]+aVals[cols-1][rows-i-1]*aVals[cols-1][rows-i-1];  
		}
		for (int j=1; j< cols; j++)
   		{	  
			a00Vals[j][0] = a00Vals[j-1][0]+aVals[j][0]*aVals[j][0];
			a10Vals[j][0] = a10Vals[j-1][0]+aVals[j][rows-1]*aVals[j][rows-1];
			a01Vals[j][0] = a01Vals[j-1][0]+aVals[cols-j-1][0]*aVals[cols-j-1][0];
			a11Vals[j][0] = a11Vals[j-1][0]+aVals[cols-j-1][rows-1]*aVals[cols-j-1][rows-1];  
			for (int i = 1; i < rows; i++)
			{
				a00Vals[j][i] = a00Vals[j][i-1]+a00Vals[j-1][i]-a00Vals[j-1][i-1]+aVals[j][i]*aVals[j][i];
				a10Vals[j][i] = a10Vals[j][i-1]+a10Vals[j-1][i]-a10Vals[j-1][i-1]+aVals[j][rows-i-1]*aVals[j][rows-i-1];
				a01Vals[j][i] = a01Vals[j][i-1]+a01Vals[j-1][i]-a01Vals[j-1][i-1]+aVals[cols-j-1][i]*aVals[cols-j-1][i];
				a11Vals[j][i] = a11Vals[j][i-1]+a11Vals[j-1][i]-a11Vals[j-1][i-1]+aVals[cols-j-1][rows-i-1]*aVals[cols-j-1][rows-i-1];  
			} 
		}
		double[][] b = new double[cols+w][rows+w];
		int rowCount = 0; int colCount = 0;

		  

		fft2D.realForwardFull(aVals);
		double[][] out = multFFT(tFFT, aVals, cols+w, 2*(rows+w));	// multiply fft_temp and fft_a
		fft2D.complexInverse(out, true);
		colCount = 0;
		for (int c=0; c<cols+w; c++)
		{
			rowCount = 0;
			for (int r=0; r<2*(rows+w); r+=2)
			{
				b[colCount][rowCount] = out[c][r];
				rowCount++;
			}
			colCount++;
		}
		int[] xy = new int[2];
		double minner = Double.POSITIVE_INFINITY;
		double minner_new = 0;
		for (int r = rows-w-1; r < rows+w; r++)
		{
			for (int c = cols-w-1; c < cols+w; c++)
			{
				if (r < rows-1 && c < cols-1)
				{
					minner_new = (a11Vals[c][r] + t00Vals[c][r] - 2*b[2*(cols-1)-c][2*(rows-1)-r])/(r*c);
				}
				else if (r >= rows-1 && c < cols-1)
				{
				    minner_new = (a01Vals[c][2*(rows-1)-r] + t10Vals[c][2*(rows-1)-r] - 2*b[2*(cols-1)-c][2*(rows-1)-r])/((2*rows-r)*c);
				}
				else if (r < rows-1 && c >= cols-1)
				{
				  	minner_new = (a10Vals[2*(cols-1)-c][r] + t01Vals[2*(cols-1)-c][r] - 2*b[2*(cols-1)-c][2*(rows-1)-r])/(r*(2*cols-c));
				}
				else 
				{
				   	minner_new = (a00Vals[2*(cols-1)-c][2*(rows-1)-r] + t11Vals[2*(cols-1)-c][2*(rows-1)-r] - 2*b[2*(cols-1)-c][2*(rows-1)-r])/((2*rows-r)*(2*cols-c));
				}

				if (minner_new < minner)
				{
					xy[0] = c; xy[1] = r; 
					minner = minner_new;
				}
			}
		}
		xy[0] = xy[0] - (cols-1);
		xy[1] = xy[1] - (rows-1);
		  
		return xy;
	}


	
	
	
	  
	  
	  
	  public static double[][] multFFT(double[][] fftOne, double[][] fftTwo, int cols, int rows)
	  {
		  double[][] product = new double[cols][rows];
		  
		  for (int c = 0; c < cols; c++)
		  {
			  for (int r = 0; r < rows; r+=2)
			  {
				  product[c][r] = fftOne[c][r]*fftTwo[c][r] - fftOne[c][r+1]*fftTwo[c][r+1];
				  product[c][r+1] = fftOne[c][r]*fftTwo[c][r+1] + fftTwo[c][r]*fftOne[c][r+1];

			  }
		  }
		  return product;
	  }
	  
	  

		public static double[] meanAndStdDev(ShortProcessor image, int rows, int cols)
		{
			 double[] ms = new double[2];
			 double mean = 0;
			 double std = 0;
			 double sum = 0;
			 double sum_sq = 0;
			 for (int r=0; r<rows; r++)
			 {
				 for (int c=0; c<cols; c++)
				 {
					 sum += image.get(c+cols*r);
					 sum_sq += image.get(c+cols*r)*image.get(c+cols*r);
				 }
			 }
			 mean = sum/(rows*cols);
			 std = Math.sqrt(sum_sq/(rows*cols) - mean*mean);

			 ms[0] = mean;
			 ms[1] = std;
			 return ms;		 
		 }	
		
		public static double[] meanAndStdDev(ByteProcessor image, int rows, int cols)
		{
			 double[] ms = new double[2];
			 double mean = 0;
			 double std = 0;
			 double sum = 0;
			 double sum_sq = 0;
			 for (int r=0; r<rows; r++)
			 {
				 for (int c=0; c<cols; c++)
				 {
					 sum += image.get(c+cols*r);
					 sum_sq += image.get(c+cols*r)*image.get(c+cols*r);
				 }
			 }
			 mean = sum/(rows*cols);
			 std = Math.sqrt(sum_sq/(rows*cols) - mean*mean);

			 ms[0] = mean;
			 ms[1] = std;
			 return ms;		 
		 }	
	  

	 public static ShortProcessor applyAffine(int[] xy, ShortProcessor oldImage, int rows, int cols)
	 {
		 ShortProcessor newImage = new ShortProcessor(cols,rows);
		 
		 if (xy[0] >= 0 && xy[1] >= 0)
		 {
			 for (int r = xy[1]; r < rows; r++)
			 {
				 for (int c = xy[0]; c < cols; c++)
				 {
					 newImage.set(c + cols*r, oldImage.get(c-xy[0] + cols*(r-xy[1])));
				 }
				 
			 }
		 }
		 
		 else if (xy[0] >= 0 && xy[1] < 0)
		 {
			 for (int r = 0; r < rows + xy[1]; r++)
			 {
				 for (int c = xy[0]; c < cols; c++)
				 {
					 newImage.set(c + cols*r, oldImage.get(c-xy[0] + cols*(r-xy[1])));
					 
				 }
			 }
			 
		 }
		 
		 else if (xy[0] < 0 && xy[1] >= 0)
		 {
			 for (int r = xy[1]; r < rows; r++)
			 {
				 for (int c = 0; c < cols+xy[0]; c++)
				 {
					 newImage.set(c + cols*r, oldImage.get(c-xy[0] + cols*(r-xy[1])));
				 }
			 }
			 
		 }
		 
		 else
		 {
			 for (int r = 0; r < rows + xy[1]; r++)
			 {
				 for (int c = 0; c < cols + xy[0]; c++)
				 {
					 newImage.set(c + cols*r, oldImage.get(c-xy[0] + cols*(r-xy[1])));
				 }
			 }
			 
		 }
		 
		 
		 return newImage;
	 }
	 
	 public static ByteProcessor applyAffine(int[] xy, ByteProcessor oldImage, int rows, int cols)
	 {
		 ByteProcessor newImage = new ByteProcessor(cols,rows);
		 
		 if (xy[0] >= 0 && xy[1] >= 0)
		 {
			 for (int r = xy[1]; r < rows; r++)
			 {
				 for (int c = xy[0]; c < cols; c++)
				 {
					 newImage.set(c + cols*r, oldImage.get(c-xy[0] + cols*(r-xy[1])));
				 }
				 
			 }
		 }
		 
		 else if (xy[0] >= 0 && xy[1] < 0)
		 {
			 for (int r = 0; r < rows + xy[1]; r++)
			 {
				 for (int c = xy[0]; c < cols; c++)
				 {
					 newImage.set(c + cols*r, oldImage.get(c-xy[0] + cols*(r-xy[1])));
					 
				 }
			 }
			 
		 }
		 
		 else if (xy[0] < 0 && xy[1] >= 0)
		 {
			 for (int r = xy[1]; r < rows; r++)
			 {
				 for (int c = 0; c < cols+xy[0]; c++)
				 {
					 newImage.set(c + cols*r, oldImage.get(c-xy[0] + cols*(r-xy[1])));
				 }
			 }
			 
		 }
		 
		 else
		 {
			 for (int r = 0; r < rows + xy[1]; r++)
			 {
				 for (int c = 0; c < cols + xy[0]; c++)
				 {
					 newImage.set(c + cols*r, oldImage.get(c-xy[0] + cols*(r-xy[1])));
				 }
			 }
			 
		 }
		 
		 
		 return newImage;
	 }

	 

	 public static int[] moveByOne2(int[] xy, ShortProcessor a, ShortProcessor t)
	 {
		 int rows = a.getHeight();
		 int cols = a.getWidth();
		 int[] newXY = new int[2];
		 int cx = Math.round(rows/250);
		 int cy = Math.round(cols/150);
		 int[] rxb = new int[cx];
		 int[] rxe = new int[cx];
		 int[] ryb = new int[cy];
		 int[] rye = new int[cy];
		 for (int i = 0; i < cx; i++)
		 {
			 rxb[i] = (int) Math.round(i*rows/((double)cx)+1);
			 rxe[i] = (int) Math.round((i+1)*rows/((double)cx));			 
		 }
		 for (int i = 0; i < cy; i++)
		 {
			 ryb[i] = (int) Math.round(i*cols/((double) cy)+1);
			 rye[i] = (int) Math.round((i+1)*cols/((double)cy));
		 }
		 double[][] z = new double[3][3];
		 
		 for (int py=0; py<cy; py++)
		 {
			 for (int px=0; px<cx; px++)
			 {
				 for (int i=-1; i<=1; i++)
				 {
					 for (int j=-1; j<=1; j++)
					 {
						 int x = i + xy[1];
						 int y = j + xy[0];
						 int boolx = (x<0)? 1 : 0; 
						 int boolx2 = (x>0)? 1 : 0;
						 int booly = (y<0)? 1 : 0;
						 int booly2 = (y>0)? 1 : 0;
						 int Lx = Math.max(rxb[px], -boolx*x+1);
						 int Rx = Math.min(rxe[px], rows-boolx2*x);
						 int Ly = Math.max(ryb[py], -booly*y+1);
						 int Ry = Math.min(rye[py], cols-booly2*y);
						 for (int r = Lx-1; r<Rx; r++)
						 {
							 for (int c = Ly-1; c < Ry; c++)
							 {
								 z[i+1][j+1] = z[i+1][j+1] + (a.get(c + cols*r) - t.get(c+y + cols*(r+x)))*(a.get(c + cols*r) - t.get(c+y + cols*(r+x)));
							 }
						 }
						  
					 }
				 }
			 }
		 }
		 for (int i = -1; i <= 1; i++)
		 {
			 for (int j = -1; j <= 1; j++)
			 {
				 z[i+1][j+1] = z[i+1][j+1]/((rows - Math.abs(i + xy[1]))*(cols - Math.abs(j + xy[0])));
			 }
		 }
		 
		 int xmin = 0;
		 int ymin = 0;

		 double minner = Double.POSITIVE_INFINITY;
		 for (int i = -1; i <= 1; i++)
		 {
			 for (int j = -1; j <= 1; j++)
			 {
				 if (z[i+1][j+1] < minner)
				 {
					 minner = z[i+1][j+1];
					 xmin = i;
					 ymin = j;
				 }
			 }
		 }
		 newXY[1] = xmin + xy[1];
		 newXY[0] = ymin + xy[0];
		
		 return newXY;
	 }
	 
	 public static int[] moveByOne2(int[] xy, ByteProcessor a, ByteProcessor t)
	 {
		 int rows = a.getHeight();
		 int cols = a.getWidth();
		 int[] newXY = new int[2];
		 int cx = Math.round(rows/250);
		 int cy = Math.round(cols/150);
		 int[] rxb = new int[cx];
		 int[] rxe = new int[cx];
		 int[] ryb = new int[cy];
		 int[] rye = new int[cy];
		 for (int i = 0; i < cx; i++)
		 {
			 rxb[i] = (int) Math.round(i*rows/((double)cx)+1);
			 rxe[i] = (int) Math.round((i+1)*rows/((double)cx));			 
		 }
		 for (int i = 0; i < cy; i++)
		 {
			 ryb[i] = (int) Math.round(i*cols/((double) cy)+1);
			 rye[i] = (int) Math.round((i+1)*cols/((double)cy));
		 }
		 double[][] z = new double[3][3];
		 
		 for (int py=0; py<cy; py++)
		 {
			 for (int px=0; px<cx; px++)
			 {
				 for (int i=-1; i<=1; i++)
				 {
					 for (int j=-1; j<=1; j++)
					 {
						 int x = i + xy[1];
						 int y = j + xy[0];
						 int boolx = (x<0)? 1 : 0; 
						 int boolx2 = (x>0)? 1 : 0;
						 int booly = (y<0)? 1 : 0;
						 int booly2 = (y>0)? 1 : 0;
						 int Lx = Math.max(rxb[px], -boolx*x+1);
						 int Rx = Math.min(rxe[px], rows-boolx2*x);
						 int Ly = Math.max(ryb[py], -booly*y+1);
						 int Ry = Math.min(rye[py], cols-booly2*y);
						 for (int r = Lx-1; r<Rx; r++)
						 {
							 for (int c = Ly-1; c < Ry; c++)
							 {
								 z[i+1][j+1] = z[i+1][j+1] + (a.get(c + cols*r) - t.get(c+y + cols*(r+x)))*(a.get(c + cols*r) - t.get(c+y + cols*(r+x)));
							 }
						 }
						  
					 }
				 }
			 }
		 }
		 for (int i = -1; i <= 1; i++)
		 {
			 for (int j = -1; j <= 1; j++)
			 {
				 z[i+1][j+1] = z[i+1][j+1]/((rows - Math.abs(i + xy[1]))*(cols - Math.abs(j + xy[0])));
			 }
		 }
		 
		 int xmin = 0;
		 int ymin = 0;

		 double minner = Double.POSITIVE_INFINITY;
		 for (int i = -1; i <= 1; i++)
		 {
			 for (int j = -1; j <= 1; j++)
			 {
				 if (z[i+1][j+1] < minner)
				 {
					 minner = z[i+1][j+1];
					 xmin = i;
					 ymin = j;
				 }
			 }
		 }
		 newXY[1] = xmin + xy[1];
		 newXY[0] = ymin + xy[0];
		
		 return newXY;
	 }




	  

}
