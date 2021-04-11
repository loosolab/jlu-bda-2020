import { Component, OnInit } from '@angular/core';
import { Router } from '@angular/router';

@Component({
  selector: 'app-home',
  templateUrl: './home.component.html',
  styleUrls: ['./home.component.scss']
})
export class HomeComponent implements OnInit {

  public graphlist: any
  
  constructor(
    private router: Router,
    ) { }

  

  ngOnInit(): void {
  }
  
  toResults(){
    this.router.navigate(["/graph_home"])
  }
  
   

}
