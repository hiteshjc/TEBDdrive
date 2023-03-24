function makeSwapgate(sites, i1,i2)
    #A swap gate that will swap site with index i1 with site with index i2
    s1 = sites[i1] #defining site with index i1
    s2 = sites[i2] #defining site with index i2
    a = ITensor(dag(s1),prime(s2))
    b = ITensor(dag(s2),prime(s1))
    for j in eachval(s1)
        a[dag(s1)[j],prime(s2)[j]] = 1.0
        b[dag(s2)[j],prime(s1)[j]] = 1.0
    end
    gate = a*b
    
    return gate
end

function swapSequence(sites, ind_1, ind_2)
    #= Producing swap gate(s) to make sure that site with index ind_1 will be right next to site with index ind_2
    The swap will be conducted on every sites, one by one
    e.g to produce swap index between 1 and 5, this function will produce a sequence of swap gate as below:
    swap (1,2), swap(2,3), swap(3,4).
    When swap(1,2) is applied, the gate will change site 1 -> site 2, and site 2 -> site 1. So now, the original
    site 1' is placed in site 2. Next, the gate will apply swap(2,3), and using the same process, the site 2 - site 3
    and site 3 -> site 2, while the original site 1', which previously placed on site 2, is now placed on site 3
    =#
    
    gate = ITensor[] #To save all the swap gates
    num_swap = abs(ind_1 - ind_2) - 1 #To fidn the number of swap needed to bring site(ind_1) and site(ind_2) next to each other
    if num_swap != 0
        if ind_2 < ind_1
            ind_1 = ind_2
        end
        for i in ind_1:(num_swap+ind_1) - 1
            U = makeSwapgate(sites, i, i+1)
            push!(gate, U)
        end
    end
    return gate
end
