package main

import (
	"fmt"
	"math"
	"math/rand"
)

type Particle struct {
	Position []float64
	Velocity []float64
	PBest    []float64
}

func NewParticle(dim int, bounds []float64) *Particle {
	p := &Particle{
		Position: make([]float64, dim),
		Velocity: make([]float64, dim),
		PBest:    make([]float64, dim),
	}

	for i := 0; i < dim; i++ {
		p.Position[i] = rand.Float64()*(bounds[1]-bounds[0]) + bounds[0]
		p.Velocity[i] = 0.0
		p.PBest[i] = p.Position[i]
	}

	return p
}

func (p *Particle) UpdatePosition() {
	for i := 0; i < len(p.Position); i++ {
		p.Position[i] += p.Velocity[i]
	}
}

func (p *Particle) UpdateVelocity(c1, c2 float64, gbest []float64) {
	r1 := rand.Float64()
	r2 := rand.Float64()

	for i := 0; i < len(p.Velocity); i++ {
		p.Velocity[i] += c1*r1*(p.PBest[i]-p.Position[i]) + c2*r2*(gbest[i]-p.Position[i])
	}
}

func (p *Particle) UpdatePBest() {
	if fitnessFunction(p.Position) < fitnessFunction(p.PBest) {
		copy(p.Position, p.PBest)
	}
}

func fitnessFunction(x []float64) float64 {
	return math.Pow(x[0], 2) + math.Pow(x[1], 2)
}

func PSO(numParticles, dim int, bounds []float64, c1, c2 float64, maxIter int) []float64 {
	particles := make([]*Particle, numParticles)
	gbest := make([]float64, dim)

	for i := 0; i < numParticles; i++ {
		particles[i] = NewParticle(dim, bounds)
	}
	copy(particles[0].Position, gbest)

	for iter := 0; iter < maxIter; iter++ {
		for _, particle := range particles {
			particle.UpdateVelocity(c1, c2, gbest)
			particle.UpdatePosition()
			particle.UpdatePBest()

			if fitnessFunction(particle.PBest) < fitnessFunction(gbest) {
				copy(particle.PBest, gbest)
			}
		}
	}

	return gbest
}

func main() {
	numParticles := 20
	dim := 2
	bounds := []float64{-10, 10}
	c1 := 2.0
	c2 := 2.0
	maxIter := 100

	gbest := PSO(numParticles, dim, bounds, c1, c2, maxIter)

	fmt.Println("Melhor posição encontrada:", gbest)
}