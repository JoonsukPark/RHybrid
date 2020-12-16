n_sim = 1000
doors = 1:3
success = 0

for (i in 1:n_sim) {
  # 차의 위치를 고른다
  car = sample(doors, 1)

  # 염소들을 배치한다
  if(car == 1) goat = c(2,3)
  else if (car == 2) goat = c(1,3)
  else goat = c(1,2)

  # 참가자가 문을 고른다
  pick = sample(doors, 1)

  # 참가자가 고르지 않은 문 중 염소가 있는 문
  goat_not_picked = goat[goat != pick]

  # 진행자가 고르지 않은 문 중 염소가 있는 문 하나를 열어준다
  if (length(goat_not_picked) > 1) open = sample(goat_not_picked, 1)
  else open = goat_not_picked

  # 참가자는 고른 문을 바꾼다
  pick = doors[(doors != pick) & (doors != open)]

  # 바꾼 문이 차가 있는 문이면 '성공'
  if (pick == car) success = success + 1
}

# '성공'한 비율
success / n_sim
