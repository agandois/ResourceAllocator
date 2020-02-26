package main;
import java.util.Set;

import input.ProblemInstance;
import output.Printer;
import solver.CplexSolver;
import solver.SCPSolver;
import solver.UserResourceInstanceAllocation;

public class Main {
	
	/**
	 * 
	 * 
	 * Command line example:
	 * PREFERENCE_FILE:/export/home/vanhee/Téléchargements/file_of_preferences.txt MIN_NB_USER_PER_RESOURCE:2 MAX_NB_USER_PER_RESOURCE:2 RESOURCE_DUPLICATE_MODE:FILE_BASED("/export/home/vanhee/Téléchargements/all_projects.csv",0,4) RESOURCE_OWNERSHIP_MODE:FILE_BASED("/export/home/vanhee/Téléchargements/all_projects.csv",0,3) PREFERENCE_MEANING:PERSONAL_INSATISFACTION OWNER_DESIRES:AT_LEAST_ONE_INSTANCE_PER_OWNER OUTPUT_MODE:LATEX_REPORT	 *  
	 * Java parameters example:
	 * -ea -Djava.library.path=/export/home/vanhee/Documents/software/cplex/cplex_studio/cplex/bin/x86-64_linux
	 *  
	 * 
	 * Contents of input.csv:
	 * RESOURCE_TYPE;USERS;PREFERENCE
	 * USERS->USER | USER,USERS
	 * RESOURCE_TYPE: any string without special characters
	 * USER: any string without special characters
	 * PREFERENCE: any integer 
	 * 
	 * Such a line means "the users grades this resource type with the given preference".
	 * The higher preference, the more preferred/important.
	 * In the program, preferences are made relative (preferred, second preferred...).
	 * 
	 * Contents of projects_by_tutor.csv:
	 * NAME;RESOURCES*
	 * RESOURCES->RESOURCE_TYPE;RESOURCES
	 * RESOURCES->RESOURCE_TYPE
	 * RESOURCE_TYPE is a resource as defined in input.csv
	 * NAME is any string without special characters
	 * 
	 * This file indicate the set of resource providers and the resource types
	 * they provide. 
	 * 
	 * Contents of nb_proj.csv:
	 * RESOURCE;AMOUNT
	 * RESOURCE is a resource as defined in input.csv
	 * AMOUNT in a positive
	 * 
	 * This file indicates the number of resource instance per resource type.
	 *  
	 * The program computes an allocation such that minimizes the insatisfaction
	 * of the most unsatisfied users while:
	 * --having at least MIN_NB_USER_PER_RESOURCE per resource instance
	 *  (except for the "remainder" if any) (to be expanded to uneven group size)
	 * --having at most MAX_NB_USER_PER_RESOURCE per resource
	 * --consuming less resource instances than the number offered by a resource provider
	 * --minimizing the consumption per resource provider (to be expanded)
	 * --never splitting groups
	 * 
	 * The preference insatisfaction is represented by the tuple
	 * (p1, p2, ... pn) where pk means that n users were given their kth preference.
	 * The aim is to first minimize pn, then given a minimal pn, minimize p(n-1)...
	 * Two interpretations for user insatisfaction are given:
	 * PERSONAL_INSATISFACTION i.e. the degree of satisfaction of the resource type 
	 * acquired by the user
	 * COMPARATIVE_INSATISFACTION i.e. the amount of people who got preferred resource 
	 * instances
	 * 
	 * Personal insatisfaction considers the relative satisfaction of the user, regardless
	 * of the allocations of other users.
	 * Comparative insatisfaction integrates the abundance of a (denied) resource type
	 * as an (aggravating) factor on user satisfaction. Being denied an abundant 
	 * resource is considered as "more unfair".
	 * 
	 * Programmer's notice: the current implementation uses mathematical tricks for
	 * finding the optimal allocation as a "one-shot" optimization operation.
	 * This causes the program to create a tight bound on the number of maximal insatisfaction,
	 * users and projects.
	 * This issue can be fixed, using a different approach (with a very similar program
	 * to what is already available).
	 * First, find the smallest pn such that the allocation is feasible. 
	 * Then, on a new program, set as a constraint that exactly pn projects should be 
	 * allocated for rank n and seek to minimize p(n-1).
	 * Proceed like this until setting all the values.
	 *  
	 * PREFERENCE_FILE:/export/home/vanhee/Téléchargements/file_of_preferences.txt MIN_NB_USER_PER_RESOURCE:2 MAX_NB_USER_PER_RESOURCE:2 RESOURCE_DUPLICATE_MODE:FILE_BASED("/export/home/vanhee/Téléchargements/all_projects.csv",0,4) RESOURCE_OWNERSHIP_MODE:FILE_BASED("/export/home/vanhee/Téléchargements/all_projects.csv",0,3) PREFERENCE_MEANING:PERSONAL_INSATISFACTION OWNER_DESIRES:AT_LEAST_ONE_INSTANCE_PER_OWNER OUTPUT_MODE:LATEX_REPORT
	 *  
	 * @param args
	 */	
	public static void main(String args[])
	{
		ProblemInstance input = ProblemInstance.newInstance(args);
				
		Set<UserResourceInstanceAllocation> res = SCPSolver.getOptimalAllocation(input);
		
		//Printer.displayResults(args,input,res);
	}

	

}