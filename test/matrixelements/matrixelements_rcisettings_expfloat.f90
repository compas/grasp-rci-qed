program matrixelements_rcisettings_expfloat
    use grasp_rciqed_rcisettings
    implicit none

    ! Write the expfloat.toml file. This will be checked by a follow-up test.
    integer :: toml_unit
    open(newunit=toml_unit, file="expfloat.toml")
    call write_toml_expfloat(toml_unit, "fermi_a", 1.23e20_real64)
    call write_toml_expfloat(toml_unit, "fermi_a", 1.23e10_real64)
    call write_toml_expfloat(toml_unit, "fermi_a", 1.23e5_real64)
    call write_toml_expfloat(toml_unit, "fermi_a", 1.23e2_real64)
    call write_toml_expfloat(toml_unit, "fermi_a", 1.23e1_real64)
    call write_toml_expfloat(toml_unit, "fermi_a", 1.23e0_real64)
    call write_toml_expfloat(toml_unit, "fermi_a", 0._real64)
    call write_toml_expfloat(toml_unit, "fermi_a", -0._real64)
    call write_toml_expfloat(toml_unit, "fermi_a", 1._real64)
    call write_toml_expfloat(toml_unit, "fermi_a", 2._real64)
    call write_toml_expfloat(toml_unit, "fermi_a", -2._real64)
    call write_toml_expfloat(toml_unit, "fermi_a", 1.23e-1_real64)
    call write_toml_expfloat(toml_unit, "fermi_a", 1.23e-2_real64)
    call write_toml_expfloat(toml_unit, "fermi_a", 1.23e-3_real64)
    call write_toml_expfloat(toml_unit, "fermi_a", 1.23e-5_real64)
    call write_toml_expfloat(toml_unit, "fermi_a", 1.23e-10_real64)
    call write_toml_expfloat(toml_unit, "fermi_a", 1.23e-20_real64)
    close(toml_unit)

end program matrixelements_rcisettings_expfloat
