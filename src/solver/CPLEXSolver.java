package solver;
import java.util.Collection;
import java.util.HashMap;
import java.util.Map;
import java.util.Optional;
import java.util.Set;
import java.util.SortedSet;
import java.util.TreeSet;
import java.util.function.Function;
import java.util.stream.Collectors;
import ilog.concert.IloException;
import ilog.concert.IloIntExpr;
import ilog.concert.IloIntVar;
import ilog.concert.IloNumExpr;
import ilog.cplex.IloCplex;
import ilog.cplex.IloCplex.UnknownObjectException;
import ilog.cplex.IloCplexModeler;
import input.ProblemInstance;
import input.ProblemInstance.OwnerDesire;
import model.ResourceOwner;
import model.User;
import model.UserGroup;
import output.Printer;

public class Solver {
	
	public enum ConstraintsOnResourceOwners{
		NO_CONSTRAINTS,	MINIMIZE_RESOURCE_CONSUMMED_PER_OWNER;}


	public static Optional<Set<UserResourceInstanceAllocation>> 
	optimizeAccordingToMaxInsatisfaction(
			int maxInsatisfaction,
			ProblemInstance inF,
			int maxNbResourcePerOwner)
	{
	
		try {
			IloCplex cplex = new IloCplex();
			cplex.setOut(null);

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
			
			Map<UserResourceInstanceAllocation, IloIntVar> allocToVar = 
					generateVariablePerAllocation(cplex,
					allAdmissibleAllocations
					);
			
			
			Map<ResourceInstance, IloIntVar> resourceMinLiftingJokerVars = 
					allAdmissibleResourceInstances.stream().collect(
							Collectors.toMap(Function.identity(),
									x->
									{
										try {
											return cplex
													.boolVar("jokerForDiscardingTheMinAllocationConstraint("+x+")"
													);
										} catch (IloException e) {
											// TODO Auto-generated catch block
											e.printStackTrace();
											throw new Error();
										}
									}));
	
			for(User pl:users)
			{
				int worseAllocValue = 1;
				if(allAdmissibleAllocations.stream()
						.anyMatch(x-> x.getUser().equals(pl)))
				{
					UserResourceInstanceAllocation worseAlloc = 
							allAdmissibleAllocations.stream()
							.filter(x-> x.getUser().equals(pl))
							.max((x,y)->inF.getInsatisfactionFor(x)-
									inF.getInsatisfactionFor(y)).get();
					worseAllocValue = inF.getInsatisfactionFor(worseAlloc);
				}
					
				/*@Deprecated
				 * for(UserResourceInstanceAllocation a: allAdmissibleAllocations)
					if(!allAdmissibleAllocations.contains(a) && a.getUser().equals(pl))
						inF.getAllocations().put(a, worseAllocValue+1);*/
			}
			
			Map<ResourceInstance, IloIntVar> varPerResource = getVarsPerResource(
					cplex, allAdmissibleResourceInstances);
			
			Map<ResourceOwner, IloIntVar> varPerOwner = getVarsPerOwner(cplex,
					inF, allAdmissibleAllocations
					);
			
			cplex.addMinimize(
					Solver.getExpressionToMinimize(
							cplex, 
							inF,
							allAdmissibleAllocations,
							users,
							inF.getRelativeInsatisfactionFor(allAdmissibleAllocations),
							allocToVar, varPerOwner));
			
			
			
			generateAllocationConstraints(
					cplex,
					allAdmissibleAllocations,
					allAdmissibleResourceInstances,
					allocToVar,
					varPerResource,
					varPerOwner,
					resourceMinLiftingJokerVars,
					inF, maxNbResourcePerOwner
					);
			
	
			cplex.exportModel("mipex1.lp");
	
			if (! cplex.solve() )
				return Optional.empty();
			System.out.println("Solution status = " + cplex.getStatus());
			
			if(inF.isDebugPrint())
				debugPrint(cplex, allocToVar, varPerResource,varPerOwner);
			Set<UserResourceInstanceAllocation>s=
					processCplexResults(cplex, allocToVar);
			
			
	
			//s.sort((x,y)->roles.indexOf(x.role) - roles.indexOf(y.role));
			System.out.println("Minimizing the number of least happy people:");
			cplex.end();
			
			return Optional.of(s);
	
		}
		catch (IloException e) {
			System.err.println("Concert exception caught '" + e + "' caught");
			throw new Error();
		}
	}
	

	private static Map<ResourceOwner, IloIntVar> getVarsPerOwner(IloCplex cplex, 
			ProblemInstance inF,
			Set<UserResourceInstanceAllocation> allAdmissibleAllocations) throws IloException {
		Map<ResourceOwner, IloIntVar>res = new HashMap<>();
		
		Set<ResourceOwner> ro = allAdmissibleAllocations.stream()
		.map(x->inF.getOwner(x.getResource()))
		.collect(Collectors.toSet());
		
		for(ResourceOwner r: ro)
			res.put(r, cplex
				.boolVar("ResourceOwnerHasAtLeastOneResourceTaken("+r+")")
				);
		return res;
	}


	private static Map<ResourceInstance, IloIntVar> getVarsPerResource(
			IloCplexModeler cplex,
			Set<ResourceInstance> allAdmissibleResourceInstances) {
		Map<ResourceInstance, IloIntVar>res = new HashMap<>();
		for(ResourceInstance r: allAdmissibleResourceInstances)
		{
			try {
				res.put(r, cplex
						.boolVar("ActiveResourceVar("+r+")"
						));
			} catch (IloException e) {e.printStackTrace();throw new Error();}
		}
		return res;
	}


	private static void debugPrint(IloCplex cplex,
			Map<UserResourceInstanceAllocation, IloIntVar> allocToVar,
			Map<ResourceInstance, IloIntVar> varPerResource,
			Map<ResourceOwner, IloIntVar> varPerOwner) throws UnknownObjectException, IloException {
		
		System.out.println("Debug print");
		
		for(UserResourceInstanceAllocation URIAlloc: allocToVar.keySet())
			System.out.println(URIAlloc+" "+cplex.getValue(allocToVar.get(URIAlloc)));
		
		for(ResourceInstance ri: varPerResource.keySet())
			System.out.println(ri+" "+cplex.getValue(varPerResource.get(ri)));
		
		for(ResourceOwner ro:varPerOwner.keySet())
			System.out.println(ro+" "+cplex.getValue(varPerOwner.get(ro)));
		
		System.out.println("Debug print end");
	}


	private static void addConstraintsOnResourceOwners(
			IloCplex cplex, 
			SortedSet<ResourceInstance> allAdmissibleResources, 
			Set<UserResourceInstanceAllocation> allAdmissibleAllocations,
			ProblemInstance inF, 
			int maxNbResourcePerOwner, 
			Map<UserResourceInstanceAllocation, IloIntVar> varPerAlloc,
			Map<ResourceInstance,IloIntVar> varPerResource) throws IloException {

		for(ResourceOwner rp : inF.getAllResourceOwners())
		{
			IloNumExpr countResourcesAllocatedForRP = cplex.constant(0);

			for(ResourceInstance r: inF.getResourceInstancesFrom(rp)
					.stream()
					.filter(x->allAdmissibleResources.contains(x))
					.collect(Collectors.toSet()))
				countResourcesAllocatedForRP = cplex.sum(
						countResourcesAllocatedForRP, varPerResource.get(r));
			cplex.addGe(
					maxNbResourcePerOwner, 
					countResourcesAllocatedForRP, 
					"AllocateResourceFromAllocatorAtMost"+maxNbResourcePerOwner+"Times("
							+rp+","+inF.getMinNumUsersPerResource()+")");
		}
	}
	


	
	

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
	 * Just a safety check
	 * @param optimal
	 * @param input
	 */
	private static void checkAllocation(
			Set<UserResourceInstanceAllocation> optimal,
			ProblemInstance input) {
		for(UserGroup ug: input.getUserGroups())
		{
			Map<User,ResourceInstance> allocatedResourcesFor=
					optimal.stream()
					.filter(x->ug.getUsers().contains(x.getUser()))
					.collect(Collectors.toMap(x->x.getUser(),
							x->x.getResource())
							);
			assert(allocatedResourcesFor.values().stream()
					.collect(Collectors.toSet()).size()==1);
		}
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

			res= Solver.optimizeAccordingToMaxInsatisfaction(
					i,
					input,Integer.MAX_VALUE);
		}

		return res.get();
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
		
		return Solver.optimizeAccordingToMaxInsatisfaction(
				worseInsatisfactionConsidered,
				input,Integer.MAX_VALUE).get();
	}



	private static boolean matchesScore(int worseInsatisfactionToMeet, int max, ProblemInstance input) {
		Optional<Set<UserResourceInstanceAllocation>> res = Optional.empty();

		System.out.println(
				"Trying to find an allocation with a maximum "
						+ "least satisfaction of rank:"+worseInsatisfactionToMeet+
						" and a maximum number of allocations per owner of:"
						+max
				);

		res= Solver.optimizeAccordingToMaxInsatisfaction(
				worseInsatisfactionToMeet,
				input,max);

		return res.isPresent();
	}


	static void generateHardAllocationsConstraints(
			IloCplex cplex, 
			Map<UserResourceInstanceAllocation, IloIntVar> varPerAlloc,
			Set<UserResourceInstanceAllocation> hardConstraints) throws IloException {
		for(UserResourceInstanceAllocation hc: hardConstraints)
			cplex.addEq(1,
					varPerAlloc.get(hc), "HardConstraint("+hc+")");
	}



	static Map<UserResourceInstanceAllocation, IloIntVar> generateVariablePerAllocation(
			IloCplex cplex,
			Collection<UserResourceInstanceAllocation> consideredAllocations
			) throws IloException {
		Map<UserResourceInstanceAllocation, IloIntVar> varPerAlloc = new HashMap<>();
		
		for(UserResourceInstanceAllocation pl : consideredAllocations)
				varPerAlloc.put(pl,
						cplex.boolVar(
						//cplex.numVar(0, 1,
								pl.toString()));
		
		return varPerAlloc;
	}



	static Set<UserResourceInstanceAllocation> processCplexResults(
			IloCplex cplex, 
			Map<UserResourceInstanceAllocation, IloIntVar> varPerAlloc) {
		return varPerAlloc.keySet().stream()
		.filter(x->{
			try {
				return cplex.getValue(varPerAlloc.get(x))>0.01;
			} catch (IloException e) {
				e.printStackTrace();
				throw new Error();
	
			}
		})
		.collect(Collectors.toSet());
	}
	
	public static void generateAllocationConstraints(
			IloCplex cplex, 
			Set<UserResourceInstanceAllocation> allAdmissibleAllocations, 
			SortedSet<ResourceInstance> allAdmissibleResources,
			Map<UserResourceInstanceAllocation, IloIntVar> varPerAlloc,
			Map<ResourceInstance, IloIntVar> allocatedResourceVar,
			Map<ResourceOwner, IloIntVar> varPerActiveOwner, 
			Map<ResourceInstance, IloIntVar> allocationJoker,
			ProblemInstance inF,
			int maxNbResourcePerOwner
			)
					throws IloException{
		
		connectResourcesAndAllocations(cplex,
				inF, 
				allocatedResourceVar,
				varPerAlloc
				);
		
		matchAllocationsOfGroups(cplex, inF, varPerAlloc);
		
		Set<UserResourceInstanceAllocation>validAllocations = varPerAlloc.keySet();
		for(User pl: inF.getAllUsers())
		{
			IloNumExpr oneResourcePerUser = 
					cplex.constant(0);
			for(UserResourceInstanceAllocation al: inF.getResouceInstanceAllocationsFor(pl))
			{
				if(!validAllocations.contains(al))continue;
				oneResourcePerUser = 
						cplex.sum(
						oneResourcePerUser,
						varPerAlloc.get(al));
			}
			cplex.addEq(
					1,
					oneResourcePerUser,
					"EachUserIsGivenExactlyOneResource("+pl+")");
		}

		
		allocateEachResourceInstanceAtMostKTimes(cplex, 
				inF, 
				varPerAlloc, 
				allAdmissibleResources,
				validAllocations);
		
		
		forceResourcesToBeAllocatedAccordingToUserResourceAllocations(
				cplex,
				allAdmissibleAllocations, allocatedResourceVar, varPerAlloc);
		
		
		forceOwnersToBeInactiveActiveIfAllOfTheirResourceAreInactive(
				cplex, inF,
				allAdmissibleAllocations, allocatedResourceVar, varPerActiveOwner);
		
		for(ResourceInstance resource: allAdmissibleResources)
		{
			IloNumExpr countUsersPerResource = cplex.constant(0);
			for(UserResourceInstanceAllocation s: inF.getAllocationsForResource(resource)
					.stream()
					.filter(x->allAdmissibleAllocations.contains(x))
					.collect(Collectors.toSet()))
			{
				countUsersPerResource = cplex.sum(
						countUsersPerResource, varPerAlloc.get(s));
			}
			
			//if the resource is not allocated, then the minimum is zero
			countUsersPerResource = cplex.sum(
					countUsersPerResource,
					cplex.prod(-inF.getMaxNbUsersPerResource(),
							allocatedResourceVar.get(resource)));
			
			countUsersPerResource = cplex.sum(
					countUsersPerResource,
					cplex.prod(-inF.getMaxNbUsersPerResource(),
							allocationJoker.get(resource)));
					
			cplex.addLe(
					countUsersPerResource,
					0, 
					"EachNonAllocatedResourceIsAllocatedAtMost0TimesUnlessTheJokerIsUsed("+resource+")");
		}
		
		
		forceAllocatedResourcesToBeAllocatedKTimes(
				cplex,
				allAdmissibleResources, inF,allAdmissibleAllocations,varPerAlloc,allocatedResourceVar, allocationJoker);
		
		
		
		
		IloNumExpr allJokers = cplex.constant(0);
		for(IloIntVar joker: allocationJoker.values())
			allJokers = cplex.sum(allJokers, joker);
		cplex.addLe(allJokers, 1,
				"AtMostOneJoker");
		
		for(ResourceInstance resource: allAdmissibleResources.stream()
				.filter(x->allAdmissibleAllocations.contains(x))
				.collect(Collectors.toSet()))
		{
			IloNumExpr countUsersPerResource = cplex.constant(0);
			for(UserResourceInstanceAllocation s: inF.getAllocationsForResource(resource))
			{
				countUsersPerResource = cplex.sum(
						countUsersPerResource, varPerAlloc.get(s));
			}
			
			//if the resource is not allocated, then the minimum is zero
			countUsersPerResource = cplex.sum(
					countUsersPerResource,
					cplex.prod(-inF.getMinNumUsersPerResource(),
							allocatedResourceVar.get(resource)));
					
			cplex.addGe(
					0, 
					countUsersPerResource, 
					"EachAllocatedResourceIsAllocatedAtLeastKTimes("
					+resource+","+inF.getMinNumUsersPerResource()+")");
		}
		
		


		//ensure pairs
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

					IloNumExpr exprUser0PicksR =
							varPerAlloc.get(u0PicksR);

					IloNumExpr exprUser1PicksR =
							varPerAlloc.get(u1PicksR);

					cplex.addEq(exprUser0PicksR, 
							exprUser1PicksR,
							"PairedUsersMatchTheirDecisions("
									+u0+","
									+u1+","+
									resource+")");
				}
			}
		
		addConstraintsOnResourceOwners(
				cplex,
				allAdmissibleResources,
				allAdmissibleAllocations,
				inF,
				maxNbResourcePerOwner, varPerAlloc, allocatedResourceVar);
		
		
		//ensure adequate resource consumption
		/*for(String resource:resources)
		{
			IloIntVar usingRKTimes = cplex.boolVar();
			IloNumExpr total = cplex.prod(usingRKTimes, inF.numberOfUsersPerResource);
			for(Allocation a: 
				varPerAlloc.keySet().stream()
				.filter(x->x.resource.equals(resource))
				.collect(Collectors.toSet()))
			{
				total = cplex.sum(total, varPerAlloc.get(a));
			}
			cplex.addEq(total, inF.numberOfUsersPerResource, "FullAllocationOfResource("+resource+")");
		}*/
		
		generateHardAllocationsConstraints(cplex, varPerAlloc, inF.getHardConstraints());

		
	}


	private static void forceOwnersToBeInactiveActiveIfAllOfTheirResourceAreInactive(
			IloCplex cplex, ProblemInstance inF,
			Set<UserResourceInstanceAllocation> allAdmissibleAllocations,
			Map<ResourceInstance, IloIntVar> allocatedResourceVar, 
			Map<ResourceOwner, IloIntVar> varPerActiveOwner) throws IloException {
		Set<ResourceInstance>instances = allAdmissibleAllocations.stream()
				.map(x->x.getResource())
				.collect(Collectors.toSet());
		Map<ResourceOwner, Set<ResourceInstance>> activeInstances = 
				inF.getOwnershipForInstances(instances);
		
		for(ResourceOwner ro: activeInstances.keySet())
		{
			IloNumExpr sumAllocatedResourceInstances = cplex.constant(0);
			for(ResourceInstance ri: activeInstances.get(ro))
			{
				sumAllocatedResourceInstances = cplex.sum(
						sumAllocatedResourceInstances, allocatedResourceVar.get(ri));
			}
			cplex.addLe(varPerActiveOwner.get(ro), sumAllocatedResourceInstances);
		}
		
		
	}


	private static void forceAllocatedResourcesToBeAllocatedKTimes(
			IloCplexModeler cplex, SortedSet<ResourceInstance> allAdmissibleResources,
			ProblemInstance inF, Set<UserResourceInstanceAllocation> allAdmissibleAllocations, Map<UserResourceInstanceAllocation, IloIntVar> varPerAlloc, Map<ResourceInstance, IloIntVar> allocatedResourceVar, Map<ResourceInstance, IloIntVar> allocationJoker) throws IloException {
		for(ResourceInstance resource: allAdmissibleResources)
		{
			IloNumExpr sumAllocatedUsersForResource = cplex.constant(0);
			for(UserResourceInstanceAllocation ua: inF.getAllocationsForResource(resource).stream()
					.filter(x->allAdmissibleAllocations.contains(x))
					.collect(Collectors.toSet()))
			{
				sumAllocatedUsersForResource = cplex.sum(
						sumAllocatedUsersForResource, varPerAlloc.get(ua));
			}
			
			//if the resource is not allocated, then the minimum is zero else it is
			//the min number of user per resource
			IloNumExpr compareAgainst = 			
			cplex.prod(inF.getMinNumUsersPerResource(), allocatedResourceVar.get(resource));
			
			cplex.addGe(
					sumAllocatedUsersForResource,
					compareAgainst, 
					"EachAllocatedResourceMustBeAllocatedAtLeastKTimes("+resource+")");
		}
	}


	private static void forceResourcesToBeAllocatedAccordingToUserResourceAllocations(
			IloCplexModeler cplex, 
			Set<UserResourceInstanceAllocation> allAdmissibleAllocations,
			Map<ResourceInstance, IloIntVar> allocatedResourceVar,
			Map<UserResourceInstanceAllocation, IloIntVar> varPerAlloc) throws IloException {
		for(UserResourceInstanceAllocation a:allAdmissibleAllocations)
		{
			cplex.addGe(
					allocatedResourceVar.get(a.getResource()),
					varPerAlloc.get(a),
					"AllocationsEntailResourceAllocations("+a.getResource()+","+a+")");
		}
	}


	private static void allocateEachResourceInstanceAtMostKTimes(
			IloCplexModeler cplex, ProblemInstance inF,
			Map<UserResourceInstanceAllocation, IloIntVar> varPerAlloc,
			Set<ResourceInstance> allAdmissibleResources,
			Set<UserResourceInstanceAllocation> allAdmissibleAllocations) throws IloException {
		for(ResourceInstance resource: allAdmissibleResources)
		{
			IloNumExpr nbOfUsersPerResource = cplex.constant(0);
			for(UserResourceInstanceAllocation ua: inF.getAllocationsForResource(resource)
					.stream()
					.filter(x->allAdmissibleAllocations.contains(x))
					.collect(Collectors.toSet()))
			{
				
				nbOfUsersPerResource = cplex.sum(
						nbOfUsersPerResource, varPerAlloc.get(ua));
			}
			cplex.addGe(
					inF.getMaxNbUsersPerResource(), 
					nbOfUsersPerResource, 
					"EachResourceIsAllocatedAtMostKTimes("+resource+","+
					inF.getMaxNbUsersPerResource()+")");
		}
	}


	private static void matchAllocationsOfGroups(IloCplex cplex, ProblemInstance inF,
			Map<UserResourceInstanceAllocation, IloIntVar> varPerAlloc
			) throws IloException {
		
		Map<UserGroup,Map<ResourceInstance,IloIntVar>> varPerResourcePerGroup = 
				new HashMap<>();

		for(UserGroup ug: inF.getUserGroups())
		{
			
			varPerResourcePerGroup.put(ug, new HashMap<>());
			for(ResourceInstance r: 
				varPerAlloc.keySet()
				.stream()
				.filter(x->ug.getUsers().contains(x.getUser()))
				.map(x->x.getResource())
				.collect(Collectors.toSet())
					)
			{
				varPerResourcePerGroup.get(ug).put(r, cplex.boolVar(
						"isResourceAllocatedToGroup("+
								r+","+ug+")"));
			}
			for(UserResourceInstanceAllocation ua:
				varPerAlloc.keySet().stream()
				.filter(x->ug.getUsers().contains(x.getUser()))
				.collect(Collectors.toSet()))
			{
				cplex.addEq(varPerResourcePerGroup.get(ug).get(ua.getResource()),
						varPerAlloc.get(ua),
						"allMemberOfUserGroupsMustBeAllocatedTheSameResource("+ua+","+ug+")"
						);
			}
		}
	}


	private static void connectResourcesAndAllocations(IloCplex cplex, ProblemInstance inF,
			Map<ResourceInstance, IloIntVar> allocatedResourceVar, 
			Map<UserResourceInstanceAllocation, IloIntVar> varPerAlloc) throws IloException 
	{
		for(ResourceInstance r: allocatedResourceVar.keySet())
		{
			IloNumExpr allAllocations = cplex.constant(0);
			
			Set<UserResourceInstanceAllocation> allocationsForResource = 
					varPerAlloc.keySet()
					.stream()
					.filter(x->x.getResource().equals(r))
					.collect(Collectors.toSet());
					
			for(UserResourceInstanceAllocation ua:allocationsForResource)
				allAllocations = cplex.sum(allAllocations, varPerAlloc.get(ua));
			
			IloNumExpr resource = 
					cplex.prod(allocationsForResource.size(), allocatedResourceVar.get(r));
			
			cplex.addLe(allAllocations, resource,
					"AllocationsMakeResourceConsummed("+r+")");
		}
	}


	public static IloNumExpr getExpressionToMinimize(
			IloCplex cplex, 
			ProblemInstance inF, 
			Set<UserResourceInstanceAllocation> allowedAllocations, 
			Set<User> users,
			Map<UserResourceInstanceAllocation, Integer> prefsPerAllocation,
			Map<UserResourceInstanceAllocation, IloIntVar> varPerAlloc, 
			Map<ResourceOwner, IloIntVar> varPerOwner
			) throws IloException {
		IloNumExpr exprToOptimize =  cplex.constant(0);
		for(UserResourceInstanceAllocation a:allowedAllocations)
		{
			exprToOptimize = 
					cplex.sum(
							exprToOptimize,
							cplex.prod(
									Math.pow(users.size(), prefsPerAllocation.get(a))+1,
									varPerAlloc.get(a)));
		}
		IloNumExpr[] weightedAllocations = new IloIntExpr[prefsPerAllocation.size()];

		int i = 0;
		for(UserResourceInstanceAllocation a: varPerAlloc.keySet())
		{
			weightedAllocations[i] = cplex.prod(
					varPerAlloc.get(a),
					prefsPerAllocation.get(a));
			i++;
		}

		if(inF.getOwnerAllocationPreferences()
				.equals(OwnerDesire.AT_LEAST_ONE_INSTANCE_PER_OWNER))
		{
			exprToOptimize = cplex.prod(varPerOwner.size()+1,
					exprToOptimize);
			
			IloNumExpr ownerInterest = cplex.constant(0);
			
			for(ResourceOwner ro:varPerOwner.keySet())
				ownerInterest = cplex.sum(cplex.prod(-1,varPerOwner.get(ro)), ownerInterest);
			
			exprToOptimize = cplex.sum(ownerInterest, exprToOptimize);	
		}
		return exprToOptimize;
	}


}
