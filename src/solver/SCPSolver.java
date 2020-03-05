package solver;

import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collection;
import java.util.Collections;
import java.util.HashMap;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;
import java.util.Map.Entry;
import java.util.Optional;
import java.util.Set;
import java.util.SortedSet;
import java.util.TreeSet;
import java.util.function.Function;
import java.util.stream.Collectors;

import input.ProblemInstance;
import input.ProblemInstance.OwnerDesire;
import model.ResourceOwner;
import model.User;
import model.UserGroup;
import scpsolver.constraints.LinearBiggerThanEqualsConstraint;
import scpsolver.constraints.LinearEqualsConstraint;
import scpsolver.constraints.LinearSmallerThanEqualsConstraint;
import scpsolver.lpsolver.LinearProgramSolver;
import scpsolver.lpsolver.SolverFactory;
import scpsolver.problems.LinearProgram;

public class SCPSolver implements ISolver {
	static Map<String , Integer> NameToId = new HashMap<>();
	
	public static Set<UserResourceInstanceAllocation> getOptimalAllocation(ProblemInstance input) {
	
		int maximumInsatisfaction = 
				SatisfactionMeasure.newInstance(
						findAllocationWithMinimalWorseInsatisfaction(input),
						input).getWorseAllocationValue();
		
		return getOptimalAllocationMatchingSatisfactionMeasureAndMinimizingResourceOwnerLoad(
				input,
				maximumInsatisfaction);
	
	}
	
	/**
	 * Computes the optimal allocation, using mathematical trickery for the reward function.
	 * Tries to minimize the maximum number of resource allocated for each owner 
	 * @param input
	 * @param worseInsatisfactionConsidered
	 * @return
	 */
	private static Set<UserResourceInstanceAllocation> getOptimalAllocationMatchingSatisfactionMeasureAndMinimizingResourceOwnerLoad(ProblemInstance input,
			int worseInsatisfactionConsidered) {	
	/*	int minAllocationPerOwner = 1;
		int maxAllocationPerOwner = 1;
		
		while(!matchesScore(worseInsatisfactionConsidered, maxAllocationPerOwner, input))
		{
			minAllocationPerOwner = maxAllocationPerOwner;
			maxAllocationPerOwner*=2;
		}
		
		Optional<Set<UserResourceInstanceAllocation>> res = Optional.empty();
		
		if(input.isMinimizingTheWorkloadOfTheMostLoaded())
		while(minAllocationPerOwner < maxAllocationPerOwner)
		{
			int i = (maxAllocationPerOwner + minAllocationPerOwner) / 2;
			System.out.println(
					"Trying to find an allocation with a maximum "
					+ "least satisfaction of rank:"+worseInsatisfactionConsidered+
					" and a maximum number of allocations per owner of:"
					+i
					);
			
			res= Solver.optimizeAccordingToMaxInsatisfaction(
					worseInsatisfactionConsidered,
					input,i);
			
			if(res.isPresent())
				maxAllocationPerOwner = i;
			else
				minAllocationPerOwner = i+1;
		}*/
		
		return SCPSolver.optimizeAccordingToMaxInsatisfaction(
				worseInsatisfactionConsidered,
				input,Integer.MAX_VALUE).get();
	}

	
	/**
	 * Find a feasible allocation, with a maximal level of insatisfaction of minimum value
	 * This allocation is not necessarily optimal.
	 * @param input
	 * @return
	 */
	private static Set<UserResourceInstanceAllocation> findAllocationWithMinimalWorseInsatisfaction(ProblemInstance input)
	{
		Optional<Set<UserResourceInstanceAllocation>> res = Optional.empty();

		for(int i = 1 ; i <= input.getAllUsers().size()&&!res.isPresent(); i++)
		{
			System.out.println(
					"Trying to find an allocation with a maximum "
							+ "least satisfaction of rank:"+i);

			res= SCPSolver.optimizeAccordingToMaxInsatisfaction(
					i,
					input,Integer.MAX_VALUE);
			System.out.println("res found");
		}

		return res.get();
	}
	
	public static Optional<Set<UserResourceInstanceAllocation>> 
	optimizeAccordingToMaxInsatisfaction(
			int maxInsatisfaction,
			ProblemInstance inF,
			int maxNbResourcePerOwner)
	{
		NameToId = new HashMap<>();
		Set<UserResourceInstanceAllocation> allAdmissibleAllocations = 
				inF.getAllocationsFilteredBy(maxInsatisfaction);

		SortedSet<User> users = 
				new TreeSet<>((x,y)->x.toString().compareTo(y.toString()));
		users.addAll(inF.getAllUsers());

		SortedSet<ResourceInstance> allAdmissibleResourceInstances = 
				new TreeSet<>((x,y)->x.toString().compareTo(y.toString()));

		allAdmissibleResourceInstances.addAll(
				allAdmissibleAllocations
				.stream().map(x->x.getResource())
				.collect(Collectors.toSet()));
			
			
		// variables booléennes
		Map<UserResourceInstanceAllocation, String> allocToVar = 
				generateVariablePerAllocation(
				allAdmissibleAllocations
				);
			
		// variables booléennes
		Map<ResourceInstance, String> resourceMinLiftingJokerVars = 
			allAdmissibleResourceInstances.stream().collect(
				Collectors.toMap(Function.identity(), x-> {
					addVar("jokerForDiscardingTheMinAllocationConstraint("+x+")");
					return "jokerForDiscardingTheMinAllocationConstraint("+x+")";
				})
			);
		
		
		// variables booléennes
		Map<ResourceInstance, String> varPerResource = getVarsPerResource(allAdmissibleResourceInstances);
		
		// variables booléennes
		Map<ResourceOwner, String> varPerOwner = getVarsPerOwner(inF, allAdmissibleAllocations);
		
		Map<UserGroup,Map<ResourceInstance,String>> varPerResourcePerGroup = new HashMap<>();
		
		for(UserGroup ug: inF.getUserGroups())
		{
			
			varPerResourcePerGroup.put(ug, new HashMap<>());
			for(ResourceInstance r: 
				allocToVar.keySet()
				.stream()
				.filter(x->ug.getUsers().contains(x.getUser()))
				.map(x->x.getResource())
				.collect(Collectors.toSet())
					)
			{
				addVar("isResourceAllocatedToGroup("+r+","+ug+")");
				varPerResourcePerGroup.get(ug).put(r, "isResourceAllocatedToGroup("+r+","+ug+")");
			}
		}
		
		
		double[] weights = SCPSolver.getExpressionToMinimize( 
				inF,
				allAdmissibleAllocations,
				users,
				inF.getRelativeInsatisfactionFor(allAdmissibleAllocations),
				allocToVar, varPerOwner);
		System.out.println(NameToId);
		for (Map.Entry<String,  Integer> entry : NameToId.entrySet()) {
			System.out.println(entry.getKey() + " " + entry.getValue());
		}
		
		LinearProgram lp = new LinearProgram(weights);
		for (int i = 0; i < weights.length; i++) {
			lp.setInteger(i);
			lp.setBinary(i);	
		}
		
		lp.setMinProblem(true);
		
		generateAllocationConstraints(
				lp,
				allAdmissibleAllocations,
				allAdmissibleResourceInstances,
				allocToVar,
				varPerResource,
				varPerOwner,
				resourceMinLiftingJokerVars,
				inF, varPerResourcePerGroup, maxNbResourcePerOwner
				);
		
		String program = lp.convertToCPLEX().toString();
		
		Map<String, Integer> sortedMap = sortByValue(NameToId);
		
		for (Map.Entry<String, Integer> entry : sortedMap.entrySet()) {
			String varName = "_" + entry.getKey();
			int varIndex = entry.getValue();
			String plVarName = "x" + varIndex;
			program = program.replace(plVarName, varName);
		}
		
		
		try {
			FileWriter file = new FileWriter("output_scpsolver.lp");
			file.write(program);
			file.close();
		} catch (IOException e) {
			e.printStackTrace();
		}
		
		//System.out.println(program);
		
		LinearProgramSolver solver  = SolverFactory.newDefault();
		System.out.println(solver);
		double[] sol = solver.solve(lp);
		System.out.println(sol);
		/*StringBuilder sb = new StringBuilder();
		for (int i = 0; i < sol.length; i++) {
			sb.append(sol[i] + " ");
			if (i % 50 == 0)
				sb.append("\n");
		}
		
		System.out.println(sb);*/
		
		Set<UserResourceInstanceAllocation> s = processPLResults(sol, allocToVar);
		
		

		//s.sort((x,y)->roles.indexOf(x.role) - roles.indexOf(y.role));
		System.out.println("Minimizing the number of least happy people:");
		
		return Optional.of(s);
	}

	static Set<UserResourceInstanceAllocation> processPLResults(
			double[] sol, 
			Map<UserResourceInstanceAllocation, String> varPerAlloc) {
		return varPerAlloc.keySet().stream()
		.filter(x->{
				return sol[NameToId.get(varPerAlloc.get(x))]>0.01;
		})
		.collect(Collectors.toSet());
	}
	
	
	static Map<UserResourceInstanceAllocation, String> generateVariablePerAllocation(
			Collection<UserResourceInstanceAllocation> consideredAllocations
			) {
		Map<UserResourceInstanceAllocation, String> varPerAlloc = new HashMap<>();
		
		
		for(UserResourceInstanceAllocation pl : consideredAllocations) {
			varPerAlloc.put(pl, pl.toString());
			addVar(pl.toString());
		}
		
		return varPerAlloc;
	}
	
	public static void generateAllocationConstraints(
			LinearProgram lp, 
			Set<UserResourceInstanceAllocation> allAdmissibleAllocations, 
			SortedSet<ResourceInstance> allAdmissibleResources,
			Map<UserResourceInstanceAllocation, String> varPerAlloc,
			Map<ResourceInstance, String> allocatedResourceVar,
			Map<ResourceOwner, String> varPerActiveOwner, 
			Map<ResourceInstance, String> allocationJoker,
			ProblemInstance inF,
			Map<UserGroup,Map<ResourceInstance,String>> varPerResourcePerGroup,
			int maxNbResourcePerOwner
			) {
		
		connectResourcesAndAllocations(lp,
				inF, 
				allocatedResourceVar,
				varPerAlloc
				);

		matchAllocationsOfGroups(lp, inF, varPerAlloc, varPerResourcePerGroup);
		
		Set<UserResourceInstanceAllocation>validAllocations = varPerAlloc.keySet();
		
		for(User pl: inF.getAllUsers())
		{
			double weights[] = new double[NameToId.size()];
			for(UserResourceInstanceAllocation al: inF.getResouceInstanceAllocationsFor(pl))
			{
				if(!validAllocations.contains(al))continue;
				int id  = NameToId.get(varPerAlloc.get(al));
				weights[id] = 1.0;
			}
			
			lp.addConstraint(new LinearEqualsConstraint(weights, 1.0,"EachUserIsGivenExactlyOneResource("+pl+")"));
			
		}
		
		
		allocateEachResourceInstanceAtMostKTimes(lp, 
				inF, 
				varPerAlloc, 
				allAdmissibleResources,
				validAllocations);
		
		forceResourcesToBeAllocatedAccordingToUserResourceAllocations(
				lp,
				allAdmissibleAllocations, allocatedResourceVar, varPerAlloc);
		
		double weights[] = new double[NameToId.size()];
		for(String joker: allocationJoker.values()) {
			int id  = NameToId.get(joker);
			weights[id] = 1.0;
			}
		lp.addConstraint(new LinearSmallerThanEqualsConstraint(weights,1.0,"AtMostOneJoker"));

	
		for(ResourceInstance resource: allAdmissibleResources.stream()
				.filter(x->allAdmissibleAllocations.contains(x))
				.collect(Collectors.toSet()))
		{
			weights = new double[NameToId.size()];
			for(UserResourceInstanceAllocation s: inF.getAllocationsForResource(resource))
			{
				int  id  = NameToId.get(varPerAlloc.get(s));
				weights[id] = 1.0;
			}
			
			int id  = NameToId.get(allocatedResourceVar.get(resource));
			weights[id] = inF.getMinNumUsersPerResource();
					
			lp.addConstraint(new LinearBiggerThanEqualsConstraint(weights,0.0,"EachAllocatedResourceIsAllocatedAtLeastKTimes("
					+resource+","+inF.getMinNumUsersPerResource()+")"));
		}
	
		for(UserGroup userGroup: inF.getUserGroups())
			for(ResourceInstance resource: allAdmissibleResources.stream()
					.filter(x->allAdmissibleAllocations.contains(x))
					.collect(Collectors.toSet()))
			{
				User u0 = userGroup.getUsers().iterator().next();
				UserResourceInstanceAllocation u0PicksR=
						UserResourceInstanceAllocation.newInstance(u0, resource);
				
				if(!varPerAlloc.containsKey(u0PicksR)) continue;
				
				for(User u1: userGroup.getUsers())	
				{
					UserResourceInstanceAllocation u1PicksR= 
							UserResourceInstanceAllocation.newInstance(u1, resource);

					weights = new double[NameToId.size()];
					int id = NameToId.get(varPerAlloc.get(u0PicksR));
					weights[id] = 1.0;
					id = NameToId.get(varPerAlloc.get(u1PicksR));
					weights[id] = -1.0;
					
					lp.addConstraint(new LinearEqualsConstraint(weights,0.0,"PairedUsersMatchTheirDecisions("
							+u0+","
							+u1+","+
							resource+")"));

				}
			}

		addConstraintsOnResourceOwners(
				lp,
				allAdmissibleResources,
				allAdmissibleAllocations,
				inF,
				maxNbResourcePerOwner, varPerAlloc, allocatedResourceVar);
		
		
		generateHardAllocationsConstraints(lp, varPerAlloc, inF.getHardConstraints());
	}
	
	private static Map<ResourceInstance, String> getVarsPerResource(Set<ResourceInstance> allAdmissibleResourceInstances) {
		Map<ResourceInstance, String> res = new HashMap<>();
		
		for(ResourceInstance r: allAdmissibleResourceInstances)
		{
			addVar("ActiveResourceVar("+r+")");
			res.put(r, "ActiveResourceVar("+r+")");
		}
		return res;
	}
	
	private static Map<ResourceOwner, String> getVarsPerOwner(
			ProblemInstance inF,
			Set<UserResourceInstanceAllocation> allAdmissibleAllocations) {
		Map<ResourceOwner, String> res = new HashMap<>();
		
		Set<ResourceOwner> ro = allAdmissibleAllocations.stream()
		.map(x->inF.getOwner(x.getResource()))
		.collect(Collectors.toSet());
		
		for(ResourceOwner r: ro) {
			addVar("ResourceOwnerHasAtLeastOneResourceTaken("+r+")");
			res.put(r, "ResourceOwnerHasAtLeastOneResourceTaken("+r+")");
		}
		return res;
	}
	
	public static double[] getExpressionToMinimize(
			ProblemInstance inF, 
			Set<UserResourceInstanceAllocation> allowedAllocations, 
			Set<User> users,
			Map<UserResourceInstanceAllocation, Integer> prefsPerAllocation,
			Map<UserResourceInstanceAllocation, String> varPerAlloc, 
			Map<ResourceOwner, String> varPerOwner
			) {
		
		double weights[] = new double[NameToId.size()];
		
		
		for(UserResourceInstanceAllocation a:allowedAllocations)
		{
			int id = NameToId.get(varPerAlloc.get(a));
			weights[id] = Math.pow(users.size(), prefsPerAllocation.get(a))+1;
		}
		
		// pas utilisé?
		/*
		IloNumExpr[] weightedAllocations = new IloIntExpr[prefsPerAllocation.size()];

		int i = 0;
		for(UserResourceInstanceAllocation a: varPerAlloc.keySet())
		{
			weightedAllocations[i] = cplex.prod(
					varPerAlloc.get(a),
					prefsPerAllocation.get(a));
			i++;
		}*/
		
		if(inF.getOwnerAllocationPreferences()
				.equals(OwnerDesire.AT_LEAST_ONE_INSTANCE_PER_OWNER))
		{
			for (int i = 0; i < weights.length; i++) {
				weights[i] = weights[i] * (varPerOwner.size()+1.0);
			}
				
			
			for(ResourceOwner ro:varPerOwner.keySet()) {
				int id = NameToId.get(varPerOwner.get(ro));
				weights[id] = -1.0;
			}	
		}
		
		return weights;
	}
	
	private static void connectResourcesAndAllocations(LinearProgram lp ,ProblemInstance inF,
			Map<ResourceInstance, String> allocatedResourceVar, 
			Map<UserResourceInstanceAllocation, String> varPerAlloc) 
	{
		for(ResourceInstance r: allocatedResourceVar.keySet())
		{
			
			Set<UserResourceInstanceAllocation> allocationsForResource = 
					varPerAlloc.keySet()
					.stream()
					.filter(x->x.getResource().equals(r))
					.collect(Collectors.toSet());
					
			double weights[] = new double[NameToId.size()];
			for(UserResourceInstanceAllocation ua:allocationsForResource) {
				int id =  NameToId.get(varPerAlloc.get(ua));
				weights[id] = 1.0;
			}

			int id = NameToId.get(allocatedResourceVar.get(r));
			weights[id] = -allocationsForResource.size();
			
			lp.addConstraint(new LinearSmallerThanEqualsConstraint(weights, 0.0,"AllocationsMakeResourceConsummed("+r+")"));

		}
	}
	
	private static void matchAllocationsOfGroups(LinearProgram lp, ProblemInstance inF,
			Map<UserResourceInstanceAllocation, String> varPerAlloc, Map<UserGroup,Map<ResourceInstance,String>> varPerResourcePerGroup
			){
		


		for(UserGroup ug: inF.getUserGroups())
		{
			
			for(UserResourceInstanceAllocation ua:
				varPerAlloc.keySet().stream()
				.filter(x->ug.getUsers().contains(x.getUser()))
				.collect(Collectors.toSet()))
			{
				double weights[] = new double[NameToId.size()];
				int id = NameToId.get(varPerResourcePerGroup.get(ug).get(ua.getResource()));
				weights[id] = -1.0;
				id =  NameToId.get(varPerAlloc.get(ua));
				weights[id] = 1.0;
				lp.addConstraint(new LinearEqualsConstraint(weights, 0.0,"allMemberOfUserGroupsMustBeAllocatedTheSameResource("+ua+","+ug+")"));
			}
		}
	}
	
	private static void allocateEachResourceInstanceAtMostKTimes(
			LinearProgram lp, ProblemInstance inF,
			Map<UserResourceInstanceAllocation, String> varPerAlloc,
			Set<ResourceInstance> allAdmissibleResources,
			Set<UserResourceInstanceAllocation> allAdmissibleAllocations){
		for(ResourceInstance resource: allAdmissibleResources)
		{
			double weights[] = new double[NameToId.size()];
			for(UserResourceInstanceAllocation ua: inF.getAllocationsForResource(resource)
					.stream()
					.filter(x->allAdmissibleAllocations.contains(x))
					.collect(Collectors.toSet()))
			{
				int id = NameToId.get(varPerAlloc.get(ua));
				weights[id] = 1.0;
			}
			
			lp.addConstraint(new LinearSmallerThanEqualsConstraint(weights,inF.getMaxNbUsersPerResource(),"EachResourceIsAllocatedAtMostKTimes("+resource+","+inF.getMaxNbUsersPerResource()+")"));
		}
	}
	
	
	private static void forceResourcesToBeAllocatedAccordingToUserResourceAllocations(
			LinearProgram lp, 
			Set<UserResourceInstanceAllocation> allAdmissibleAllocations,
			Map<ResourceInstance, String> allocatedResourceVar,
			Map<UserResourceInstanceAllocation, String> varPerAlloc){
		for(UserResourceInstanceAllocation a:allAdmissibleAllocations)
		{
			double weights[] = new double[NameToId.size()];
			int id = NameToId.get(varPerAlloc.get(a));
			weights[id] = -1.0;
			
			id = NameToId.get(allocatedResourceVar.get(a.getResource()));
			weights[id] = 1.0;
			
			lp.addConstraint(new LinearBiggerThanEqualsConstraint(weights,0.0,"AllocationsEntailResourceAllocations("+a.getResource()+","+a+")"));
		}
	}
	
	
	private static void addConstraintsOnResourceOwners(
			LinearProgram lp, 
			SortedSet<ResourceInstance> allAdmissibleResources, 
			Set<UserResourceInstanceAllocation> allAdmissibleAllocations,
			ProblemInstance inF, 
			int maxNbResourcePerOwner, 
			Map<UserResourceInstanceAllocation, String> varPerAlloc,
			Map<ResourceInstance,String> varPerResource){

		for(ResourceOwner rp : inF.getAllResourceOwners())
		{
			double weights[] = new double[NameToId.size()];
			for(ResourceInstance r: inF.getResourceInstancesFrom(rp)
					.stream()
					.filter(x->allAdmissibleResources.contains(x))
					.collect(Collectors.toSet())) {
				int id = NameToId.get(varPerResource.get(r));
				weights[id] = 1.0;
			}
			
			lp.addConstraint(new LinearSmallerThanEqualsConstraint(weights,maxNbResourcePerOwner,"AllocateResourceFromAllocatorAtMost"+maxNbResourcePerOwner+"Times("
					+rp+","+inF.getMinNumUsersPerResource()+")"));
			
		}
	}
	
	static void generateHardAllocationsConstraints(
			LinearProgram lp, 
			Map<UserResourceInstanceAllocation, String> varPerAlloc,
			Set<UserResourceInstanceAllocation> hardConstraints){
		for(UserResourceInstanceAllocation hc: hardConstraints) {
			double weights[] = new double[NameToId.size()];
			int id = NameToId.get(varPerAlloc.get(hc));
			weights[id] = 1.0;
			lp.addConstraint(new LinearEqualsConstraint(weights, 1.0,"HardConstraint("+hc+")"));
		}
	}
	
	private static void addVar(String name) {
		if (!NameToId.containsKey(name))
			NameToId.put(name, NameToId.size());
	}
	
    public static <K, V extends Comparable<? super V>> Map<K, V> sortByValue(Map<K, V> map) {
        List<Entry<K, V>> list = new ArrayList<>(map.entrySet());
        list.sort(Entry.comparingByValue());

        Map<K, V> result = new LinkedHashMap<>();
        for (int i = list.size() - 1; i >= 0; i--) {
        	Entry<K, V> entry = list.get(i);
            result.put(entry.getKey(), entry.getValue());
        }

        return result;
    }
		
}
