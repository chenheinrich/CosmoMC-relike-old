    module KdeTests
    use KdeReionizationTanhHighz
    use KernelDensityEstimate
    implicit None
    integer, parameter, private :: sp= KIND(1.d0)

    contains

    function RunKdeTests() result(fails)
    integer fails
    fails = 0
    
    call test_kde(fails)
    end function RunKdeTests

    subroutine test_kde(fails)
    integer fails
    Type(TIniFile) :: Ini
    real(mcp) :: LogLike
    real(mcp), dimension(10) :: P
    
    type(RelikeKde) :: Kde 
    type(RelikeKdeChains):: KdeChains

    call RelikeKde_Init(this%Kde, this%KdeChains, 0.2_mcp)
    P(1) = 0.04_mcp !taulo
    P(2) = 0.02_mcp !tauhi
    call RelikeKde_OneModel(P, LogLike) 
    print *, 'LogLike = ', LogLike
        
    end subroutine test_kde

    end module KdeTests