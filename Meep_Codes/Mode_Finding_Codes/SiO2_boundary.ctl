(set-param! resolution 50) ; resolution 
(define-param sx 10) ; Cell size in x
(define-param sy 4) ; Cell size in y
(define-param sz 1) ; Cell size in z
(define n_bto 2.278) ; Estimated refractive index of BTO
(define-param w_len_um 1.55) ; light wavelength in um
(define-param lip_width_um 1.5) ; Lip width in um
(define-param lip_height_um 0.175) ; Lip height in um
(define-param bto_height_um 0.325) ; BTO layer height in um
(define-param sio2_height_um 2) ; height of sio2 in um
(define-param sio2_thickness_nm 50)
(define-param dpml 1) ; Thickness of pml
(define-param df 0.01) ; Frequency spread
(define-param pad 0.1234) ; Distance from source to edge of cell
(define-param k_interp 31)
(define-param output_field false) ; True or false value for outputting a field
(define-param k_z 0.6)
(define-param k_low 1.175)
(define-param k_high 1.5)
(define-param num_pics 10) ; Number of output pictures
(define-param freq (/ 1 w_len_um)) ; frequency of light in meep units
(define-param output_eps_only false) ; Parameter saying wheter or not you only want to output the eps file
(define sio2_thickness_um (/ sio2_thickness_nm 1000))

(set! geometry-lattice (make lattice (size sx sy sz)))

(set! geometry (list (make block (size infinity bto_height_um infinity) (material (make medium (index n_bto))) (center 0 (+ -1.5 (/ sio2_height_um 2) (/ bto_height_um 2)) 0))
                     (make block (size infinity sio2_height_um infinity) (material (make medium (index 1.44))) (center 0 -1.5 0))
                     (make block (size lip_width_um lip_height_um infinity) (material (make medium (index n_bto))) (center 0 (+ -1.5 (/ sio2_height_um 2) bto_height_um (/ lip_height_um 2)) 0))

		     (make block (size lip_width_um sio2_thickness_um infinity) (material (make medium (index 1.44))) (center 0 (+ -1.5 (/ sio2_height_um 2) bto_height_um lip_height_um (/ sio2_thickness_um 2)) 0))                     

))

(set! sources (list
               (make source
                 (src (make gaussian-src (frequency freq) (fwidth df)))
                 (component Hz) (direction Z) (size lip_width_um (* (+ lip_height_um bto_height_um) 0.95) 0) (center 0 (+ -1.5 (/ sio2_height_um 2) (/ (+ bto_height_um lip_height_um) 2)) (+ (* -0.5 sz) pad)))

))

(set! pml-layers (list (make pml (direction Y) (thickness dpml))))
(set! pml-layers (list (make pml (direction X) (thickness dpml))))


(if output_eps_only
  (begin
  (init-fields)
  (output-epsilon)
  )
  (begin
  (if output_field 
    (begin
    (set-param! k-point (vector3 0 0 k_z))
    (run-sources+ 300 (at-beginning output-epsilon))
    (run-until (/ 1 freq) 
      (at-every (/ 1 freq num_pics) output-dpwr)
    )
    )
    (begin
    (run-k-points 300 (interpolate k_interp (list (vector3 0 0 k_low) (vector3 0 0 k_high))))
    )
  )
  )
)

(exit)
